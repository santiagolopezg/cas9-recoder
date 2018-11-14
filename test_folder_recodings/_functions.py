# functions used by recode_1.py and indels.py

'''
To do:
- calculate distance from cutsite and use as reference for best pam;
- ??
'''



import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, IUPAC
import math


#load reference sequence
def load_ref_seq(file):
	with open(file, 'r') as f:
		data = f.read()
	if '.ape' in file: # formats input .ape file to a string seq output
		data = data.split('ORIGIN\r\n')[1]
		seq = ''
		for i in data:
			if i in ['A','T','C','G','N']:
				seq = seq+i 
	if '.txt' in file:
		seq = data
	return seq

# translate seq to protein
def translate(ref_seq):
	seq = Seq(ref_seq[ref_seq.find('ATG'):])
	transcribed_seq = seq.transcribe()
	#print transcribed_seq
	translated_refseq =  transcribed_seq.translate(to_stop=True)
	#print translated_refseq
	return translated_refseq

def mk_aa_pos(translated_refseq):
	#print translated_refseq
	aa = []
	pos = []
	c = 1
	for i in translated_refseq:
		aa.append(i)
		pos.append(c)
		c+=1
	#print aa_pos_dict
	return aa, pos
		
		

#get position of first start codon
def position_atg(refseq, u_d_stream):
	pos = refseq.find('ATG')
	#print pos
	d = 1 - pos #normalize to put ATG at 1
	#print d
	return pos, d
	
	
# returns both aa and possible codons, when given only one of them
# this function is for the codon_change.py 
def get_codon_and_aa(input, pos, d, ref_seq, u_d_stream, aminoacids, codon_table, flag):
	#input is the abbreviation of the aminoacid / letter representation
	#pos is aa position while the actual positions from get_positions will be 3*(pos-1)

	#print 'position', pos
	#print pos
	l = input
	if input in aminoacids.keys():
		l = aminoacids[input]
	codons = codon_table[l]
	#print codons
	for codon in codons:
		positions = get_positions(codon, ref_seq, u_d_stream, d, flag)
		#print pos
		#print positions
		if pos in positions:
			return codon, l, pos
	print "Couldn't find the amino acid in that position :-("
	exit()
		
# scan sequence for all sequence occurrences
def get_positions(codon, ref_seq, u_d_stream, d, flag):
	positions = []
	index = 0
	while index < len(ref_seq):
		index = ref_seq.find(codon, index)
		#print index
		if index == -1:
			break
		#print '- {0} found at {1}'.format(codon, index-u_d_stream+1)
		
		# check whether there are enough nucleotides to upstream
		# region of the codon on sequence
		
		p = index + d
		if flag == 'aa':
			p= (p-1)/3 +1
		elif flag == 'indels':
			p = p-1
		#print p
		positions.append(p)
		index += 1 
	#print positions
	return positions

# get up and downstream nucleotide sequences to see where I can anchor the guide
def check_around(ref_seq, codon, position, u_d_stream , d, flag):
	
	positions = get_positions(codon, ref_seq, u_d_stream , d, flag)
	#print positions
	if position not in positions:
		print 'Wrong sequence / position. Please re-try.'
		exit()
	else:
		#print d, position, position - d
		#print position
		if flag == 'aa':
			npoo = ((position-1)*3 +1) - d
			#print npoo
			n_pos = npoo
		if flag == 'indels':
			#position is relative to ATG -- turn it into absolute
			n_pos = position - d
			# npos should be 32
			#n_pos = position
		#print "Foodle", n_pos
		l_limit = n_pos-u_d_stream
		u_limit = n_pos+u_d_stream
		#print l_limit, u_limit
		if l_limit < 0:
			#position += l_limit 	# this so that relative position isn't 
			l_limit = 0			# distorted if position is not in the middle of region
		if u_limit > len(ref_seq):
			u_limit = len(ref_seq)
		region = ref_seq[l_limit:u_limit]
		return region, position, l_limit, u_limit
	
# returns protein sequence corresponding to region of interest, but starting
# from ATG and spanning len(region interest)-pos_1st_ATG
def p_seq_region(translated_seq, region):
	# since we're allocating 30 bp downstream of ATG, that's 10 aa
	trans = translated_seq[:(len(region)-region.find('ATG', 100))/3]
	trans_formatted = ''
	for aa in trans:
		trans_formatted += ' '+aa+' '
	return trans_formatted
	

# finds all occurrences of PAM sequences within the region of interest
def find_pam(region, pam_dict,  pam, u_d_stream, atg, d, ref_seq, l_limit, u_limit):
	
	possible_pams = pam_dict['NGG']
	print 'possible PAMs: {0}\n'.format(possible_pams)
	pam_positions = {} #dictionary that maps PAM -> positions where it's found
	for pams in possible_pams: #goes through all possible pams and tries to find them in refseq
		index = 0
		pam_positions[pams]=[]
		while index < len(region):
			index = region.find(pams, index)
			if index == -1:
				break
			#print '- {0} found at {1}'.format(pams, index-u_d_stream+1)
			
			#this needs revising - something wrong with pam location is calculated
			
			p = index + d + l_limit
			#poo = (p-1)/3 + 2
			
			#print pams, index, atg, p, d, l_limit, u_limit

			#print poo
			pam_positions[pams].append(p)
			#pam_positions[pams].append(index-u_d_stream+1)
			index += 1
	#print pam_positions
	return pam_positions


# returns best PAM, and its position
def best_pam_pos(pam_positions, u_d_stream, pos, atg, d, ref_seq, l_limit, u_limit, flag):
	max_d = 2*u_d_stream
	for i in pam_positions.keys(): #for each pam
		for j in pam_positions[i]: # for its position (relative or absolute - depends)
			if flag == 'indels':
				dif = j-pos		# positions for indels refer to positions relative
				#dif = j
				#print dif, i
			if flag == 'aa':	# to the ATG (in nucleotides); for indels, in aminoacids
				#print i, pos
				poopoo = (pos-1)*3+1
				#print 'position of aa: ', poopoo
				if 'GG' in i: 				#AGTATGGACAAGAAGTACT
												#in this scenario, distance of cutsite to the first A in AAG is 9
					dif = j - poopoo -4
					print j, poopoo, dif, i

				if 'CC' in i:
					dif = j - poopoo +5
				#print dif
			if max_d >= int(math.fabs(dif)):
			
				
				#if dif >= 1 and 'GG' in i: # this makes sure that the NGG PAM is downstream
				
				if 'GG' in i:
					smallest_d = dif
					max_d = math.fabs(dif)   # of the n20 on the reference sequence
					best_pam = i
					#best_pam_d = j+u_d_stream
					best_pam_d = j - d - l_limit
					smallest_d = dif
					
					#print best_pam_d
				#if dif < 1 and 'CC' in i: # this makes sure of the same thing, just that
				if 'CC' in i:
					smallest_d = dif
					max_d = math.fabs(dif)		# since it's on the complement strand, it 
					best_pam = i				# looks upstream
					#best_pam_d = j+u_d_stream
					best_pam_d = j - d - l_limit
	#exit()
	#print "print d, best_pam, best_pam_d", d, best_pam, best_pam_d
	return best_pam, best_pam_d, smallest_d 

def make_guides(region, ref_seq, pam, pam_location, l_limit, atg):
	#print pam_location, l_limit, atg
	pam_location +=l_limit #translate local position to general position (on refseq)
	if 'GG' in pam:
		guide = ref_seq[pam_location-20:pam_location] # looks upstream of PAM to build guide
		#should return ACGTACACACAAGGCGGTAA
	if 'CC' in pam: 
		guide = ref_seq[pam_location+3:pam_location+23] # looks downstream of PAM 
		guide = str(Seq(guide, generic_dna).reverse_complement())
				
		#guide = guide[::-1]
		#should return TCACATAACTTAAGAAGGGA				# on complement strand (so downstream 3-5
	return guide									 # ~ upstream on 5-3 strand)

def store_g_d(guide, donor, codon, pos, dict_of_guides,  dict_of_donors,smallest_d_dict, after, smallest_d, flag):
	if flag == 'aa':
		#print pos
		pos = pos
	dict_of_guides['{0}>{1}'.format(codon, after)] = guide
	dict_of_donors['{0}>{1}'.format(codon, after)] = donor
	smallest_d_dict[guide]=math.fabs(smallest_d)
	#print dict_of_guides
	return dict_of_guides, dict_of_donors, smallest_d_dict

def save_guides_and_donors(dict_of_guides, dict_of_donors, smallest_d_dict, best_pam, j, file_name, pos, flag):
	with open('best_g+d_{0}_{1}_{2}_{3}.txt'.format(file_name, j, pos, flag), 'w') as f:
		f.write('seq/aa \tposition \tbest pam\t aminoacid dist from cut site \t\t guide \t\t\t donor\n')
		for i in dict_of_guides.keys():
			codon_pos = i
			guide = dict_of_guides[i]
			donor =  dict_of_donors[i]
			
			smallest_d = smallest_d_dict[guide]
			best_pam = best_pam
			#if 'CC' in best_pam:
			#	donor = str(Seq(donor, generic_dna).reverse_complement())
			f.write('{0} \t {1}\t {2} \t {3}\t {4}\t {5}\n'.format(codon_pos,pos, best_pam, smallest_d, guide, donor))
		

def make_donor(ref_seq, best_pam, best_pam_d, max_d, codon, position, after, pam_positions, atg, d, aminoacids, codon_table, donor, i, flag):		
	#takes the ref_seq, the codon and pam their respective positions and returns 
	# a donor sequence to make the desired change. len(donor) = 90; codon to change
	# should be right in the middle of it. Should also modify PAM to be unrecognisable.
	# very similar to look_around function (ref_seq, codon, position, u_d_stream , d, flag)
	#print d, position, position - d
	#print position

	if flag == 'aa':
		n_pos = ((position-1)*3 +1) - d
	if flag == 'indels':
		#position is relative to ATG -- turn it into absolute
		n_pos = position - d
		# npos should be 32
		#n_pos = position
	#print "Foodle", n_pos
	
	
	print position, n_pos, d
	#print 'BEST PAM DDDD', atg, d, best_pam_d, max_d, best_pam_d-2-80, best_pam_d-2+80-len(codon)
	abc = best_pam_d #notice best_pam_d is distance from start of region, not to abs start 
	#u_hom = n_pos-45
	u_hom = n_pos-100#+max_d-2-80
	#d_hom = n_pos+45-len(codon)
	d_hom = n_pos+97#+max_d-2+80-len(codon)
	
	
	#print u_hom, d_hom, n_pos
	
	yapa = 0
	if u_hom < 0:
		#position += l_limit 	# this so that relative position isn't 
		yapa = -1*u_hom
		u_hom = 0			# distorted if position is not in the middle of region
	
	print 'esta es la yapa', yapa
	if d_hom > len(ref_seq):
		d_hom = len(ref_seq)
		
	#print u_hom, d_hom, n_pos
	print u_hom, n_pos, d_hom
	u_hom_donor = ref_seq[u_hom:n_pos]
	l_hom_donor = ref_seq[n_pos+len(codon):d_hom+len(codon)+yapa]
	
	#print 'up and down homology: ', u_hom_donor, l_hom_donor
	

	#best_codon = mutations_to_seq(codon, n_pos, donor, ref_seq, after, aminoacids, codon_table, pam_positions, wtd = 'change aa')
	
	best_codon = after
	
	print best_codon
	if len(donor) == 0:
		donor = u_hom_donor+best_codon+l_hom_donor
		#print position, u_hom, best_codon
	else:
		#print position, u_hom, best_codon
		
		#exit()
		
		#print u_hom_donor, donor.find(u_hom_donor[60:70])
		
		donor, best_codon = list(donor), list(best_codon)
		
		#print donor[n_pos-u_hom+3*i], donor[n_pos-u_hom+3*i+1], donor[n_pos-u_hom+3*i+2], best_codon[0], best_codon[1], best_codon[2]

		donor[n_pos-u_hom+3*i], donor[n_pos-u_hom+3*i+1], donor[n_pos-u_hom+3*i+2] = best_codon[0], best_codon[1], best_codon[2]
		donor = ''.join(str(e) for e in donor)
		print donor
		
	
	#print codon, best_codon, after, n_pos, best_pam_d
	
	#codon_pam_dict['{0}_{1}'.format(codon, after)] = [n_pos, best_codon, best_pam, best_pam_d]
	
	#print codon_pam_dict
	#print position
	return donor, position
		
	
def change_pam(best_pam, best_pam_d, atg, d, position, codon_table, aminoacids, pam_positions, ref_seq, donor, wtd='change pam'):

	#change the pam, if it's != from codon changed by donor
	#print best_pam, best_pam_d, donor.find(best_pam, best_pam_d)
	
	
	#print 'oh man, gotta smash the pam'
	print 'PAM needs to be changed - will attempt to'
	
	
	if 'GG' in best_pam:
		print 'GG tho'

	elif 'CC' in best_pam:
		print 'CC tho'

	wtd = 'change pam'
	p = (position-1)*3 +1 - atg - d
	position = (position-1)*3 +1
	# if the pam is CCN, that means the cut will be on the comp. strand;
	# so what has to be done is modify comp strand codons and check them 

	#print position, position/3, position%3
	#print pam

	if p%3 == 0: #PAM is a codon
		for j in codon_table.keys():
			if best_pam in codon_table[j]:
				possible_codons = codon_table[j]
				#print possible_codons
		changed_pam = mutations_to_seq(best_pam, best_pam_d, donor, ref_seq, possible_codons, aminoacids, codon_table, pam_positions, wtd = 'change pam')
		#print changed_pam
		print 'changed pam to:', changed_pam
		donor, changed_pam = list(donor), list(changed_pam)
		donor[best_pam_d], donor[best_pam_d+1], donor[best_pam_d+2] = best_pam[0], best_pam[1], best_pam[2]
		donor = ''.join(str(e) for e in donor)
		#donor = donor.replace(best_pam, changed_pam)


	if p%3 != 0: #PAM is in between 2 codons	
		#figure out between which two codons it is
		# could be like XNG|GXX (position%3 = 1) or XXN|GGX (position%3 = 2)
	
		if p%3==1:
			#XNG|GXX -> preferably change the first G (1st codon, last nucleotide)
			#print 'meeeeeee'

			p_codon = donor.find(best_pam, best_pam_d)-1
			n_codon = donor.find(best_pam, best_pam_d)+2
			prev_codon = donor[p_codon:p_codon+3]
			next_codon = donor[n_codon:n_codon+3]
			
			#print 'fooooo', prev_codon, next_codon
			for j in codon_table.keys():
				if prev_codon in codon_table[j]:
					possible_codons = codon_table[j]
			for c in possible_codons:
				a = [c[0] for c in possible_codons]
			try:
				a = a.remove('G')
			except:
				print 'Not only Gs'				
			
			if a == None:
				for k in codon_table.keys():
					if next_codon in codon_table[k]:
						possible_codons = codon_table[k]
				#print 'fooooo', next_codon
				
				
				changed_pam = mutations_to_seq(next_codon, p_codon, donor, ref_seq, possible_codons, aminoacids, codon_table, pam_positions, wtd = 'change pam')
				print 'changed pam to:', changed_pam
				donor, changed_pam = list(donor), list(changed_pam)
				donor[best_pam_d], donor[best_pam_d+1], donor[best_pam_d+2] = changed_pam[0], changed_pam[1], changed_pam[2]
				donor = ''.join(str(e) for e in donor)
				#donor = donor.replace(prev_codon, changed_pam)

			else:
				changed_pam = mutations_to_seq(prev_codon, n_codon, donor,  ref_seq, possible_codons, aminoacids, codon_table, pam_positions, wtd = 'change pam')
				print 'changed pam to:', changed_pam
				donor, changed_pam = list(donor), list(changed_pam)
				donor[best_pam_d], donor[best_pam_d+1], donor[best_pam_d+2] = changed_pam[0], changed_pam[1], changed_pam[2]
				donor = ''.join(str(e) for e in donor)
				#donor = donor.replace(next_codon, changed_pam)				
			
		elif p%3==2:
		
			print "PAM straight in middle of codon"
			#XXN|GGX -> GAA|GGG -> preferably change the 1st G in NGG (2nd codon, 1st nucleotide)
			#TAA GAA GGG AGC
			#print best_pam, donor.find(best_pam, best_pam_d)
			p_codon = donor.find(best_pam, best_pam_d)-2
			n_codon = donor.find(best_pam, best_pam_d)+1
			prev_codon = donor[p_codon:p_codon+3]
			next_codon = donor[n_codon:n_codon+3]
			#print 'fooooo', prev_codon, next_codon
			for j in codon_table.keys():
				if next_codon in codon_table[j]:
					possible_codons = codon_table[j]
			for c in possible_codons:
				a = [c[0] for c in possible_codons]
			try:
				a = a.remove('G')
			except:
				print 'Not only Gs'				
			
			if a == None:
				for k in codon_table.keys():
					if prev_codon in codon_table[k]:
						possible_codons = codon_table[k]
				#print 'fooooo', prev_codon
				
				changed_pam = mutations_to_seq(prev_codon, p_codon, donor, ref_seq, possible_codons, aminoacids, codon_table, pam_positions, wtd = 'change pam')
				donor, changed_pam = list(donor), list(changed_pam)
				donor[best_pam_d], donor[best_pam_d+1], donor[best_pam_d+2] = changed_pam[0], changed_pam[1], changed_pam[2]
				donor = ''.join(str(e) for e in donor)
				#donor = donor.replace(prev_codon, changed_pam)

			else:
				changed_pam = mutations_to_seq(next_codon, n_codon, donor,  ref_seq, possible_codons, aminoacids, codon_table, pam_positions, wtd = 'change pam')
				donor, changed_pam = list(donor), list(changed_pam)
				donor[best_pam_d], donor[best_pam_d+1], donor[best_pam_d+2] = changed_pam[0], changed_pam[1], changed_pam[2]
				donor = ''.join(str(e) for e in donor)
				#donor = donor.replace(next_codon, changed_pam)

			
				#print changed_pam

	#print position, n_pos
	return donor

			
def mutations_to_seq(codon, n_pos, donor, ref_seq, after, aminoacids, codon_table, pam_positions, wtd = 'change pam'):
	#takes the aminoacid we want to mutate and what we want to change it to, and
	# and returns the codon for the mutated aa.
	
	if type(after) == list:
		if len(after) == 1:
			print 'Only codon corresponding to that aminoacid, cannot change identified pam ({0}) :-('.format(after[0])
			return after[0]
		after_codons = after
		
	elif len(after) == 1: #gave me an aminoacid letter
		after_codons = codon_table[after]
	elif after in aminoacids.keys():
		after_letter = aminoacids[after]
		after_codons = codon_table[after_letter]
	#print codon, after_codons
	
	#print 'pololololo', donor.find(codon, n_pos)-3, donor.find(codon, n_pos), donor.find(codon, n_pos)+3
	
	if wtd == 'change pam':
		if 'CC' in codon:
			#print 'cc', codon, donor.find(codon, n_pos), n_pos
			prev_codon = donor[donor.find(codon, n_pos)-3:donor.find(codon, n_pos)]
			next_codon = donor[donor.find(codon, n_pos)+3:donor.find(codon, n_pos)+6]
	
		if 'GG' in codon:
			#print 'gg', donor.find(codon, n_pos)
			prev_codon = donor[donor.find(codon, n_pos)-3:donor.find(codon, n_pos)]
			next_codon = donor[donor.find(codon, n_pos)+3:donor.find(codon, n_pos)+6]
	
	else:
		
		prev_codon = ref_seq[ref_seq.find(codon, n_pos)-3:ref_seq.find(codon)]
		next_codon = ref_seq[ref_seq.find(codon, n_pos)+3:ref_seq.find(codon)+6]
		
	print 'upst and downst codons: ',prev_codon, next_codon
	
	valid_codons = []
	gs = len(codon)
	dg = 0
	cs = len(codon)
	dc = 0

	best_codon = codon
	for codon in after_codons:
		#print codon
		around = prev_codon+codon+next_codon
		#print around 
		for pam in pam_positions.keys():
			if pam in around:
				break
			if codon not in valid_codons and codon not in pam_positions.keys():
				valid_codons.append(codon)
		#print valid_codons
		
		g_s = around.count('G')
		d_g = g_s - gs
		c_s = around.count('C')
		d_c = c_s - gs
		
		#print d_g, d_c
		if d_g < dg or d_c < dc:
			dg = d_g
			dc = d_c
			best_codon = codon
	#print best_codon
	return best_codon
	

	
###########################

#print functions

def print_codon(pos, codon, region, u_d_stream, atg):
	dist = atg-int(math.fabs(pos))
	return'-'*(dist)+'#'*len(codon)+'-'*(len(region)-dist-len(codon))

def print_aa(pos, codon, region, u_d_stream, atg, d, l_limit, ref_seq): #have to fix this
	#dist = atg+int(math.fabs(pos))-1
	#dist = u_d_stream - int(math.fabs(pos))
	pos = (pos-1)*3 +1
	pos = pos -(1 - ref_seq.find('ATG'))-l_limit
	#return'-'*(region.find(codon))+'#'*len(codon)+'-'*(len(region)-len(codon)-(region.find(codon)))
 	return'-'*(pos)+'#'*len(codon)+'-'*(len(region)-pos-len(codon))
 	#len(codon)-(region.find(codon)))
 
def print_start(atg, region, flag, d):
	if flag == 'aa':
		if region.find('ATG') != -1:
			return '-'*(region.find('ATG'))+'1'+'-'*(len(region)-region.find('ATG')-1)+'\n\n'
		else: 
			return '-'*len(region)
	if flag == 'indels':
		return '-'*(atg)+'1'+'-'*(len(region)-atg)+'\n\n'
	
def print_pseq_region(atg, p_seq_region, region, flag):
	if flag == 'aa':
		return '-'*(region.find('ATG'))+p_seq_region
	if flag == 'indels':
		return '-'*(atg-1)+p_seq_region
	
def print_pam_loc(best_pam_d, region):
	return '-'*(best_pam_d)+'*'+'-'*(len(region)-best_pam_d-1)
	#return '-'*(10)+'*'+'-'*(len(region)-best_pam_d)
	
	
def print_guide_loc(region, guide):
	return '-'*(region.find(guide))+guide+'-'*(len(region)-len(guide)-region.find(guide))



	
	

