''' 
change codons that encode specific aminoacids from specific coordinates

* Now, idea is to:
	-> change aa in position x to Alanine, for all aas
	-> change aa in positions x and x+1 to alanine, for all aas
	-> change aa in positions x, x+1 and x+2 to alanine for all aas

'''

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import math
from _functions import *
import random

homology_arms = int(input('Enter desired length of homology arms (l <= 100; recommended: 66): '))
if homology_arms > 100:
	print 'l must be <= 100. Please re-try.'
	exit()

file = raw_input('Enter name of reference sequence (.txt format): ')

if '.txt' not in file:
	print 'Please enter a reference sequence in .txt format, flanked by 100 bp upstream and downstream sequences.'
	exit()
#file = 'helicase.txt'

mode = int(input('Enter tiling mode (aminoacids tiled per donor - 1 / 2 / 3): '))

u_d_stream = 100
file_name = file.split('.')[0]
after = raw_input('Enter aminoacid used for tiling (default is A): ')

pam_dict = {'NGG': ['AGG', 'TGG', 'CGG', 'GGG', 'CCA', 'CCT', 'CCG', 'CCC']}
aminoacids = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
   			  'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
 		      'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
  		      'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'STOP': '*'}
  		      
codon_table = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['TTT', 'TTC'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'I': ['ATT', 'ATC', 'ATA'],
    'H': ['CAT', 'CAC'],
    'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'M': ['ATG'],
    'N': ['AAT', 'AAC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    '*': ['TAA', 'TAG', 'TGA'],
}



flag = 'aa'
pam = 'NGG'
donor = ''

######################################################


ref_seq = load_ref_seq(file)
atg, d = position_atg(ref_seq, u_d_stream)

translated_refseq = translate(ref_seq)
aa, positions = mk_aa_pos(translated_refseq)

#print aa, positions

guides = []
donors = []
events = []
pms = []
dst=[]

#print aa
for p in xrange(len(aa)-(mode-1)): #for all of the aminoacids we want to change:
	#print after, p

	guide = ''
	donor = ''
	c = p
	codons_before = []
	codons_after = []
	possible_codons_after = codon_table[after][:]
	pam = 'NGG'

	pams = {}
	
	for q in xrange(mode):
		aminoacid = aa[c+q]
		position = positions[c+q]
		
		
		codon, l, position = get_codon_and_aa(aminoacid, position, d, ref_seq, u_d_stream, aminoacids, codon_table, flag)

		codons_before.append([codon, l, position])
		
		print ('Corresponding codon is {0}.'.format(codon))
		# translate the nucleotide sequence to protein


		# retrieve sequence to look into for the guide
		region, position, l_limit, u_limit = check_around(ref_seq, codon, position, u_d_stream, d, flag)

		# find PAM positions 
		pam_positions = find_pam(region, pam_dict, pam, u_d_stream, atg, d, ref_seq, l_limit, u_limit)

		#find best pam and its position
		best_pam, best_pam_d, max_d = best_pam_pos(pam_positions, u_d_stream, position, atg, d, ref_seq, l_limit, u_limit, flag = 'aa')

		guide = make_guides(region, ref_seq, best_pam, best_pam_d, l_limit, atg)
		
		
		pams['{0}-{1}'.format(best_pam, q)] = [math.fabs(max_d), best_pam_d, guide]

		#print aminoacid, position
		
		#print q, c
		
		#print best_pam, best_pam_d, math.fabs(max_d)
				
		#c += 1
	
	#print pams
	
	max_d = 10000
	for pam_ in pams.keys():
		if pams[pam_][0] < max_d:
			best_pam = pam_.split('-')[0]
			index = int(pam_.split('-')[1])
			best_pam_d = pams[pam_][1]
			guide = pams[pam_][2]
			
			max_d = int(pams[pam_][0])
	
	#print best_pam, max_d, guide, best_pam_d, index, codons_before[index], codons_before
	
	#now we know which is the best guide, best pam, its position and etc, we just have
	# to build 1 donor which incorporates changes to 1 / 2 / 3 aminoacids > A
	
	
	print 'best pam is', best_pam
	if 'GG' in best_pam:
		print codons_before, len(codons_before)
		for i in xrange(mode):
		
		#for q in codons_before:
			
			#print 'qqqqq', q
			after_c = random.choice(possible_codons_after)
			print 'changing codon to ', after_c, codons_before[i][2]
			donor, position = make_donor(ref_seq, best_pam, best_pam_d, max_d+3, codons_before[i][0], codons_before[i][2], after_c, pam_positions, atg, d, aminoacids, codon_table, donor, i, flag)
			possible_codons_after.remove(after_c)
			
			print donor
			
			
			#exit()
			#print possible_codons_after
			#print donor, position
		
		#exit()
	if 'CC' in best_pam:
		for i in xrange(mode):
			after_c = random.choice(possible_codons_after)
			#print after
			donor, position = make_donor(ref_seq, best_pam, best_pam_d, max_d-5, codons_before[i][0], codons_before[i][2], after_c, pam_positions, atg, d, aminoacids, codon_table, donor, i, flag)
			possible_codons_after.remove(after_c)
		#donor = str(Seq(donor, generic_dna).reverse_complement())
		#modify other positions without building a new guide
	
	
	print guide, donor
	
	pms.append(best_pam)
	
	if donor.find(best_pam, best_pam_d) != best_pam_d:
		events.append('pam changed by mutation')
	
	if donor.find(best_pam, best_pam_d) == best_pam_d:
		
		print 'shieeeeet', donor.find(best_pam, best_pam_d), best_pam_d
		
		n_pos = ((position-1)*3 +1) - d
		u_hom = n_pos-100
		d_hom = n_pos+97
		yapa = 0
		if u_hom < 0:
			#position += l_limit 	# this so that relative position isn't 
			yapa = -1*u_hom
			u_hom = 0			# distorted if position is not in the middle of region

		#print 'esta es la yapa', yapa
		if d_hom > len(ref_seq):
			d_hom = len(ref_seq)
		
		c_pos = n_pos-u_hom
					 
		#print guide in donor
		print 'npos', n_pos
		#print 'donor find guide', donor.find(guide)+20
		print c_pos
		print best_pam_d # position on donor
		print c_pos - best_pam_d
		
		res = (c_pos - best_pam_d)%3
		
		ra = 3-res
		
		ppos = best_pam_d
		#need to find codon position 
		
		#exit()
		
		if 'GG' in best_pam:
		
			if res == 2:
			
				c_pos_s = ppos +5
				c_pos_l = ppos +2
				print donor[c_pos_l:c_pos_s]
				
			
			if res == 1:
			
				c_pos_s = ppos +4
				c_pos_l = ppos +1
				print donor[c_pos_l:c_pos_s]
			
			if res == 0:
			
				c_pos_s = ppos + 3
				c_pos_l = ppos
				print donor[c_pos_l:c_pos_s]
				
				
				
			#	print c_pos_l, c_pos_s
			#	c = donor[c_pos_l:c_pos_s]
			#	print c
			
			if res == 0 and best_pam == 'GGG':
				print guide in donor
				
				#exit()
				
			
			g = 'yes'
			
				#guide_pos = donor.find(guide)+21
			while g == 'yes' and (donor[best_pam_d:best_pam_d+3] in pam_dict['NGG']):
				#print 'ala moana'
				cod = donor[c_pos_l:c_pos_s]
				
				print cod
			
				for codons in codon_table.values():
					if cod == after_c:
						print 'hamlet', cod, after_c, codons
						break
					#print 'foogar'
					elif cod in codons and len(codons)>1:
					#	print 'agar'
						print cod, codons
						sample = codons[:]
						print sample
						sample.remove(cod)
						print sample
						for s in sample:
							print 'trying codon {0}'.format(s)
							donor, s = list(donor), list(s)
							donor[c_pos_l], donor[c_pos_l+1], donor[c_pos_l+2] = s[0], s[1], s[2]
							donor = ''.join(str(e) for e in donor)
						
							print donor.find(best_pam, best_pam_d) == best_pam_d
							if (donor[best_pam_d:best_pam_d+3] not in pam_dict['NGG']):
								break
							if guide not in donor:
								g = 'no'
								break
								
							print s
						if (donor[best_pam_d:best_pam_d+3] not in pam_dict['NGG']):
							break
						if guide not in donor:
							g = 'no'
							break
							
				print donor
							
							
			#	if res == 0 and best_pam == 'GGG':
					#exit()

				c_pos_s -= 3
				c_pos_l -= 3
				
			print donor
			print donor.find(best_pam, best_pam_d) == best_pam_d
			if donor.find(best_pam, best_pam_d) != best_pam_d:
				events.append('changed pam')
			if guide not in donor:
				events.append('changed guide')
			print guide in donor
				
		#exit()
		
		
		if 'CC' in best_pam:
			c_pos = n_pos-u_hom
					
		
			ppos = best_pam_d
					 
			
			#print guide in donor
			print n_pos
			print donor.find(guide)+20
			print c_pos
			print best_pam_d # position on donor
			print c_pos - best_pam_d
		
			res = (c_pos - best_pam_d)%3
		
			ra = 3-res
		
			#exit()
			if res == 2:
			
				c_pos_s = ppos +2
				c_pos_l = ppos -1
				print donor[c_pos_l:c_pos_s]
				
				#exit()
				
				
			
			if res == 1:
			
				c_pos_s = ppos +1
				c_pos_l = ppos -2
				print donor[c_pos_l:c_pos_s]
			
				
				
			if res == 0:
			
				c_pos_s = ppos + 3
				c_pos_l = ppos 
				print donor[c_pos_l:c_pos_s]
			
				#exit()
				
			print res
			
			
			
			#	print c_pos_l, c_pos_s
			#	c = donor[c_pos_l:c_pos_s]
			#	print c
			
				
			g = 'yes'
				#guide_pos = donor.find(guide)+21
			while g == 'yes' and (donor[best_pam_d:best_pam_d+3] in pam_dict['NGG']):
				cod = donor[c_pos_l:c_pos_s]
				print cod
			
				for codons in codon_table.values():
				
					if cod == after_c:
						print 'hamlet', cod, after_c, codons
						break
					elif cod in codons and len(codons)>1:
						print cod, codons
						sample = codons[:]
						print sample
						sample.remove(cod)
						print sample
						for s in sample:
							print 'trying codon {0}'.format(s)
							donor, s = list(donor), list(s)
							donor[c_pos_l], donor[c_pos_l+1], donor[c_pos_l+2] = s[0], s[1], s[2]
							donor = ''.join(str(e) for e in donor)
						
							print donor.find(best_pam, best_pam_d) == best_pam_d
							
							if (donor[best_pam_d:best_pam_d+3] not in pam_dict['NGG']):
								break
							if guide not in donor:
								g = 'no'
							 	break 

							
							print s
					
						if (donor[best_pam_d:best_pam_d+3] not in pam_dict['NGG']):
							break
						if guide not in donor: 
							g = 'no'
							break 		
					
					
				c_pos_s += 3
				c_pos_l += 3
				
			print donor
			print donor.find(best_pam, best_pam_d) == best_pam_d
			if donor.find(best_pam, best_pam_d) != best_pam_d:
				events.append('changed pam')
			if guide not in donor:
				events.append('changed guide')
			print guide in donor
			#guide_pos = donor.find(guide)-2
				
		print donor
		#exit()
	#print 'length of donor:', len(donor)
	if len(donor) != 200:
		print 'Donor misconstructed. Please report bug.'

	guides.append(guide)
	donors.append(donor)
	dst.append(math.fabs(max_d))
	
	#exit()


print len(events), len(donors)
with open('best_g+d_{0}_{1}{2}.txt'.format(file_name, mode, after), 'w') as f:	
	f.write('aminoacid change(s) \t position(s) \t PAM \t dist from cutsite to change pos \t guide\t donor\n\n')
	for i in xrange(len(guides)):
		#print i
		aminoa = ''
		p = ''
		for j in xrange(mode):
			#print j, i+j
			aminoa += aa[i+j]
		
		for k in xrange(mode):
			p += str(positions[i+k])+','
		guide = guides[i]
		pam = pms[i]
		d = dst[i]
		#event = events[i]
		donor = donors[i][int(math.fabs(u_d_stream-homology_arms)):(u_d_stream+homology_arms+3)]
		f.write('{0}>{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(aminoa, after*mode, p, pam, int(d), guide, donor))

		#f.write('\t{0}>{1}\t{2}\t{3}\t{4}\t{5}\n'.format(aminoa, 'A'*mode, p, event, guide, donor))

exit()
	