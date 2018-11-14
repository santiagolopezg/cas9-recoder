# Codon changes to all possible others

from Bio.Seq import Seq
from Bio.Alphabet import generic_rna, generic_dna
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

u_d_stream = 100

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

file_name = file.split('.')[0]

dict_of_guides = {}
dict_of_donors = {}
codon_pam_dict = {}	
smallest_d_dict = {}
				
		
########################################################


ref_seq = load_ref_seq(file)
atg, d = position_atg(ref_seq, u_d_stream)

translated_refseq = translate(ref_seq)

aa, positions = mk_aa_pos(translated_refseq)

print aa, positions


pam = 'NGG'

for p in xrange(len(aa)): ## for all aminoacids in the protein
	
	guides = []
	donors = []
	changes = []
	#events = []
	pms = []
	psps = []
	

	for aas in codon_table.keys(): ## change to all aminoacids	
	
		
		position = positions[p]
	
		guide = ''
		donor = ''
		codons_before = []
		codons_after = []
		aa_bef = aa[p]
		
		if aa_bef != aas:
		
			psps.append(p)
			changes.append(aas)
		
			after = random.choice(codon_table[aas])
		
	
			print p, aa_bef, aas, after
		
		
			codon, l, position = get_codon_and_aa(aa_bef, position, d, ref_seq, u_d_stream, aminoacids, codon_table, flag)

			#codons_before.append([codon, l, position])
		
			print ('Corresponding codon is {0}.'.format(codon))
			# translate the nucleotide sequence to protein
		
		

			# retrieve sequence to look into for the guide
			region, position, l_limit, u_limit = check_around(ref_seq, codon, position, u_d_stream, d, flag)

			# find PAM positions 
			pam_positions = find_pam(region, pam_dict, pam, u_d_stream, atg, d, ref_seq, l_limit, u_limit)

			#find best pam and its position
			best_pam, best_pam_d, max_d = best_pam_pos(pam_positions, u_d_stream, position, atg, d, ref_seq, l_limit, u_limit, flag = 'aa')

			guide = make_guides(region, ref_seq, best_pam, best_pam_d, l_limit, atg)
		
			pms.append(best_pam)
		
			print 'best pam is', best_pam 
			print guide, best_pam, best_pam_d, max_d
		
			if 'GG' in best_pam:
				#print codons_before, len(codons_before)	
				donor, position = make_donor(ref_seq, best_pam, best_pam_d, max_d+3, codon, position, after, pam_positions, atg, d, aminoacids, codon_table, donor, 0, flag)
			
			if 'CC' in best_pam:
	
				donor, position = make_donor(ref_seq, best_pam, best_pam_d, max_d-5, codon, position, after, pam_positions, atg, d, aminoacids, codon_table, donor, 0, flag)

			print donor
		
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
							if cod == after:
								print 'hamlet', cod, after, codons
								break
							#print 'foogar'
							if cod in codons and len(codons)>1:
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
									
									print donor[best_pam_d:best_pam_d+3], best_pam

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
					#if donor.find(best_pam, best_pam_d) != best_pam_d:
					#	events.append('changed pam')
					#if guide not in donor:
					#	events.append('changed guide')
					#print guide in donor
			
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
							if cod == after:
								print 'hamlet', cod, after, codons
								break
							if cod in codons and len(codons)>1:
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
					
									print (donor[best_pam_d:best_pam_d+3] in pam_dict['NGG'])
						
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
					print (donor[best_pam_d:best_pam_d+3] in pam_dict['NGG'])
				#	if donor.find(best_pam, best_pam_d) != best_pam_d:
						#events.append('changed pam')
				#	if guide not in donor:
				#		events.append('changed guide')
					print guide in donor
					#guide_pos = donor.find(guide)-2
			
				print donor
				#exit()
			
			#print 'length of donor:', len(donor)
			if len(donor) != 200:
			
				print 'Donor misconstructed. Please report bug.'
				exit()

			guides.append(guide)
			donors.append(donor)
		
	
	with open('best_g+d_{0}_{1}_{2}_aa.txt'.format(file_name, aa[p], positions[p]), 'w') as f:	
		f.write('aa \tposition \tbest pam\t aminoacid dist from cut site \t\t guide \t\t\t donor\n')
		
		#print len(guides), len(aa), len(changes), changes, changes[0]
		a = aa[p]
		for i in xrange(len(guides)):
			#print i
			change = changes[i]
			if a != change:
				p = psps[i]+1
				guide = guides[i]
				pam = pms[i]
				#event = events[i]
				#donor = donors[i][14:146]
				donor = donors[i][int(math.fabs(u_d_stream-homology_arms)):(u_d_stream+homology_arms+3)]
				#print donor
				#print len(donor)
				f.write('{0}>{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(a, change, p, pam, int(math.fabs(max_d)), guide, donor))

	#exit()	
		
		# 66 -> 14:146
		# 146 = 