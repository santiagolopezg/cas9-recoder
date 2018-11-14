'''
twofold validation of sequences:
1. align the donor with the reference sequence
2. translate it and align it with reference (also translated);
'''

import glob, os
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna, generic_dna
import math
#from _functions import *
import random
import zipfile

#load reference sequence


ref_seq = raw_input('Enter name of reference sequence (i.e. .txt file used for previous steps): ')
whattodo = raw_input('Enter type of sequences to validate {r (recoding) // t (tiling)}: ')
rs = raw_input('Enter restriction site sequence {eg. AGGCCT}: ')

name = ref_seq.split('.')[0]
ref_seq = open(ref_seq, 'r').read()

	
good = 0
bad = 0
total = 0
refs = []
posss = []
pams = []
dists = []
guides = []
donors = []

#print ref_seq
#print name
#open file with donor(s)
if whattodo == 'r':

	for file in glob.glob('best_g+d_{0}_*_*_aa_filtered.txt'.format(name)):
		#print file
		with open(file, 'r') as f:
			data = f.read().split('\n')
			data = [e for e in data if len(e) > 0]
		
		for l in data:
			l = l.split('\t')
			#print l
		
			ref, pos, pam, d, guide, donor = l[0], l[1], l[2], l[3], l[4], l[5]
		
			p = ref_seq.find(donor[:20])
			p_rseq = ref_seq[p:p+len(donor)]
		
			if 'CC' in pam:
				guide = str(Seq(guide).reverse_complement())
				print_guide = '###'+guide
			else:
				print_guide = guide+'###'
			
			p_guide = ref_seq.find(guide)
		

			## first check
		
			print 'nucleotide alignment of reference and donor\n'
			print '{0}, position {1}'.format(ref, pos)
			if 'CC' in pam:
				print ' '*(p_guide-p-3)+print_guide
			else:
				print ' '*(p_guide-p)+print_guide
			print p_rseq
			print donor
			print ''
			#print p
		

			#print ref_seq.find('ATG', 100)
		
			translated_ref = Seq(ref_seq[ref_seq.find('ATG',100):], generic_dna).transcribe().translate(to_stop=True)
		
			#now, insert donor into refseq
		
			translated_ref_with_donor = ref_seq[:p]+donor+ref_seq[p+len(donor):]
				
			print 'aminoacid alignment of reference and donor\n'

			try:
				translated_ref_p = translated_ref[int(pos)-22:int(pos)+22]
				assert len(translated_ref_p) > 0
				translated_ref_with_donor_p = Seq(translated_ref_with_donor[translated_ref_with_donor.find('ATG',100):], generic_dna).transcribe().translate(to_stop=True)
			
				translated_ref_with_donor_pee = translated_ref_with_donor_p[int(pos)-22:int(pos)+22]
				t_r_donor = Seq(translated_ref_with_donor[translated_ref_with_donor.find('ATG',100):], generic_dna).transcribe().translate(to_stop=True)
			
				print ' '*21+'#'
			except: # this happens when mutating the M in position 1, because it screws the whole downstream translation and has
				translated_ref_p = translated_ref[:44]		# no upstream coding regions (obviously)
				translated_ref_with_donor_p = Seq(translated_ref_with_donor[100:], generic_dna).transcribe().translate(to_stop=True)
				translated_ref_with_donor_pee = translated_ref_with_donor_p[:44]
				t_r_donor = translated_ref_with_donor_p

				print ' '*(int(pos)-1)+'#'

			## second check
		
			print translated_ref_p
			print translated_ref_with_donor_pee
			print ''

		
			print '-'*len(donor)+'\n'

			## third check
		
		
			translated_ref, t_r_donor = list(translated_ref), list(t_r_donor)
		
		
			#print translated_ref, t_r_donor
		
			#exit()
		
			for i in xrange(len(t_r_donor)):
				a, b = translated_ref[i], t_r_donor[i]
				#print a, b
				if a != b:
					#print a, b
					#exit()
					if b == ref.split('>')[1]:
						good +=1
					else:
						print ref.split('>')[1], 'Error. Recoding affects other aminoacids. Stopping run.'
						bad +=1
						exit()
					
					total +=1
			
	print 'good: {0}; bad: {1}; total: {2}'.format(good, bad, total)

elif whattodo == 't':
	for file in glob.glob('best_g+d_{0}_*_filtered_for_pos.txt'.format(name)):
		#print file
		with open(file, 'r') as f:
			data = f.read().split('\n')
			data = [e for e in data if len(e) > 0]
		
		for l in data:
			l = l.split('\t')
			#print l
		
			ref, pos, pam, d, guide, donor = l[0], l[1], l[2], l[3], l[4], l[5]
		
			poso = pos.split(',')[0]
			dig = len(pos.split(',')) -1 
			#print dig 
			#exit()
			p = ref_seq.find(donor[:13])
			#print p, 'pee'
			p_rseq = ref_seq[p:p+len(donor)]
		
			if 'CC' in pam:
				guide = str(Seq(guide).reverse_complement())
				print_guide = '###'+guide
			else:
				print_guide = guide+'###'
			
			p_guide = ref_seq.find(guide)
		

			## first check

			print 'nucleotide alignment of reference and donor\n'
			print '{0}, position {1}'.format(ref, pos[:len(pos)-1])
			if 'CC' in pam:
				print ' '*(p_guide-p-3)+print_guide
			else:
				print ' '*(p_guide-p)+print_guide
			print p_rseq
			print donor
			print ''
			#print p
		

			#print ref_seq.find('ATG', 100)
		
			translated_ref = Seq(ref_seq[ref_seq.find('ATG',100):], generic_dna).transcribe().translate(to_stop=True)
		
			#now, insert donor into refseq
		
			translated_ref_with_donor = ref_seq[:p]+donor+ref_seq[p+len(donor):]
				
			print 'aminoacid alignment of reference and donor\n'

			try:
				translated_ref_p = translated_ref[int(poso)-22:int(poso)+22]
				assert len(translated_ref_p) > 0
				translated_ref_with_donor_p = Seq(translated_ref_with_donor[translated_ref_with_donor.find('ATG',100):], generic_dna).transcribe().translate(to_stop=True)
			
				translated_ref_with_donor_pee = translated_ref_with_donor_p[int(poso)-22:int(poso)+22]
				t_r_donor = Seq(translated_ref_with_donor[translated_ref_with_donor.find('ATG',100):], generic_dna).transcribe().translate(to_stop=True)
			
				print ' '*21+('#'*dig)
			except: # this happens when mutating the M in position 1, because it screws the whole downstream translation and has
				translated_ref_p = translated_ref[:44]		# no upstream coding regions (obviously)
				translated_ref_with_donor_p = Seq(translated_ref_with_donor[100:], generic_dna).transcribe().translate(to_stop=True)
				translated_ref_with_donor_pee = translated_ref_with_donor_p[:44]
				t_r_donor = translated_ref_with_donor_p

				print ' '*(int(poso)-1)+'#'

			## second check
		
			print translated_ref_p
			print translated_ref_with_donor_pee
			print ''

		
			print '-'*len(donor)+'\n'

			## third check
		
		
			translated_ref, t_r_donor = list(translated_ref), list(t_r_donor)
		
		
		
			#print translated_ref, t_r_donor
		
			#exit()
		
			for i in xrange(len(t_r_donor)):
				for j in xrange(dig):
					#print i+j-1, len(translated_ref), len(t_r_donor)
					try:
						a, b = translated_ref[i+j-1], t_r_donor[i+j-1]
					except:
						break
					#print a, b
				#print a, b
					if a != b:
						#print a, b
						#exit()
						potato = ref.split('>')[1]
						#print b, potato[j-1]
						if b == ref.split('>')[1][j-1]:
							good +=1
						else:
							print ref.split('>')[1], 'Error. Recoding affects other aminoacids. Stopping run.'
							bad +=1
							exit()
					
						total +=1
			
	print 'good: {0}; bad: {1}; total: {2}'.format(good, bad, total)

		
if bad == 0:

	if whattodo == 'r':

		for file in glob.glob('best_g+d_{0}_*_*_aa_filtered.txt'.format(name)):
			#print file
			with open(file, 'r') as f:
				data = f.read().split('\n')
				data = [e for e in data if len(e) > 0]
				
			nn = file.split('.txt')[0]+'_rs.txt'
						
			with open(nn, 'w') as g:
			
				for line in data:
					l = line.split('\t')
					#print l
		
					ref, pos, pam, d, guide, donor = l[0], l[1], l[2], l[3], l[4], l[5]
										
					##troubleshooting 
					
					#rs = 'ACNNNNNCTCC' # 5 Ns
					#guide = 'POTATO'
					#donor = 'AAAAAAACNNNNNCTCCAAAAAA'
					
					rs_c = str(Seq(rs).reverse_complement())
					
					if 'N' in rs:
						numb_n = 0
						for c in rs:
							if c=='N':
								numb_n +=1
								
						rs_1 = rs.split('N')[0]
						rs_2 = rs.split('N')[-1]
						
						rs_c1 = rs_c.split('N')[0]
						rs_c2 = rs_c.split('N')[-1]
					
						#print rs_1, rs_2, rs_c1, rs_c2, numb_n
						
						if guide.find(rs_1)!=-1:
							if guide.find(rs_2)== guide.find(rs_1)+numb_n+len(rs_1):
								print 'RS found'
						elif donor.find(rs_1)!=-1:
							if donor.find(rs_2)== donor.find(rs_1)+numb_n+len(rs_1):
								print 'RS found'
						elif guide.find(rs_c1)!=-1:
							if guide.find(rs_c2)== guide.find(rs_c1)+numb_n+len(rs_c1):
								print 'RS found'
						elif donor.find(rs_c1)!=-1:
							if donor.find(rs_c2)== donor.find(rs_c1)+numb_n+len(rs_c1):
								print 'RS found'
						else:
							g.write(str(line)+'\n')
						
					else:			
						if rs not in guide and rs not in donor and rs_c not in guide and rs_c not in donor:
							g.write(str(line)+'\n')
		
		filez = []
		for f in glob.glob('best_g+d_{0}_*_*_aa_filtered_rs.txt'.format(name)):  
			if os.stat(f).st_size == 0:
		   		os.remove(f)
			else:
	   			filez.append(f)
		
		if len(filez)> 0:
			with zipfile.ZipFile('{0}_g+d_recodings_filtered_rs.zip'.format(name), 'w') as zip:
				for file in filez:
					zip.write(file)
					os.remove(file)

			
		with zipfile.ZipFile('{0}_g+d_recodings.zip'.format(name), 'w') as zip:
			for f in glob.glob('best_g+d_{0}_*_*_aa_filtered.txt'.format(name)):  
				if os.stat(f).st_size == 0:
		   			os.remove(f)
		   		else:
					zip.write(f)
			
		for file in glob.glob('best_g+d_{0}_*_*_aa_filtered.txt'.format(name)):
			os.remove(file)	
			
			
	elif whattodo == 't':
	
		for file in glob.glob('best_g+d_{0}_*_filtered_for_pos.txt'.format(name)):

			#print file
			with open(file, 'r') as f:
				data = f.read().split('\n')
				data = [e for e in data if len(e) > 0]
				
			nn = file.split('.txt')[0]+'_rs.txt'
						
			with open(nn, 'w') as g:
			
				for line in data:
					l = line.split('\t')
					#print l
		
					ref, pos, pam, d, guide, donor = l[0], l[1], l[2], l[3], l[4], l[5]
										
					##troubleshooting 
					
					#rs = 'ACNNNNNCTCC' # 5 Ns
					#guide = 'POTATO'
					#donor = 'AAAAAAACNNNNNCTCCAAAAAA'
					
					rs_c = str(Seq(rs).reverse_complement())
					
					if 'N' in rs:
						numb_n = 0
						for c in rs:
							if c=='N':
								numb_n +=1
								
						rs_1 = rs.split('N')[0]
						rs_2 = rs.split('N')[-1]
						
						rs_c1 = rs_c.split('N')[0]
						rs_c2 = rs_c.split('N')[-1]
					
						#print rs_1, rs_2, rs_c1, rs_c2, numb_n
						
						if guide.find(rs_1)!=-1:
							if guide.find(rs_2)== guide.find(rs_1)+numb_n+len(rs_1):
								print 'RS found'
						elif donor.find(rs_1)!=-1:
							if donor.find(rs_2)== donor.find(rs_1)+numb_n+len(rs_1):
								print 'RS found'
						elif guide.find(rs_c1)!=-1:
							if guide.find(rs_c2)== guide.find(rs_c1)+numb_n+len(rs_c1):
								print 'RS found'
						elif donor.find(rs_c1)!=-1:
							if donor.find(rs_c2)== donor.find(rs_c1)+numb_n+len(rs_c1):
								print 'RS found'
						else:
							g.write(str(line)+'\n')
						
					else:			
						if rs not in guide and rs not in donor and rs_c not in guide and rs_c not in donor:
							g.write(str(line)+'\n')
		
		filez = []				
			
		for f in glob.glob('best_g+d_{0}_*_filtered_for_pos_rs.txt'.format(name)):  
			if os.stat(f).st_size == 0:
		   		os.remove(f)
			else:
	   			filez.append(f)
		
		if len(filez)> 0:
			with zipfile.ZipFile('{0}_g+d_tilings_filtered_rs.zip'.format(name), 'w') as zip:
				for file in filez:
					zip.write(file)
					os.remove(file)

			
		with zipfile.ZipFile('{0}_g+d_tilings.zip'.format(name), 'w') as zip:
			for f in glob.glob('best_g+d_{0}_*_filtered_for_pos.txt'.format(name)):  
				if os.stat(f).st_size == 0:
		   			os.remove(f)
		   		else:
					zip.write(f)
			
		for file in glob.glob('best_g+d_{0}_*_filtered_for_pos.txt'.format(name)):
			os.remove(file)		
			
	print 'files zipped!'
	
