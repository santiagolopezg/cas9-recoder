# sanity check all the files and filter against 1. guides in donor sequences and 2. presence
# of PAMS in donors


import os, glob

c = 0 # count how many I have to toss out
good = 0
not_in_but_far = 0
nibtf_20 = 0
in_but_far = 0
pam_not_destr = 0
total = 0 # counts the total changes I have
nb = 0

pams_GG = ['AGG', 'TGG', 'CGG', 'GGG']
pams_CC = ['CCA', 'CCT', 'CCG', 'CCC']

whattodo = raw_input('Enter type of files to filter {r (recoding) // t (tiling)}: ')

# possible files: best_g+d_*_*_*_aa.txt // best_g+d_*_*A.txt

if whattodo == 'r':
	for file in glob.glob('best_g+d_*_*_*_aa.txt'):
		nb+= 1
		name = file.split('.')[0]
		nameee = file.split('_')[2]
		print name
		with open(file, 'r') as f:
			with open('{0}_filtered.txt'.format(name), 'w') as g:
				#g.write('seq/aa \tposition \tbest pam\t aminoacid dist from cut site \t\t guide \t\t\t donor\n')
				data = f.read().split('\n')
				data = data[1:len(data)-1]
				#print len(data), data
		
				for line in data:
					total +=1 
					line = line.split('\t')
					change, position, pam, d, guide, donor = line[0], line[1], line[2], line[3], line[4], line[5]
					#print change, position, pam, d, guide, donor
					d = int(float(d))

					if guide in donor:
						print guide, donor 
						# fuck
						if 'CC' in pam:
							reg = donor[donor.find(guide)-3:donor.find(guide)]
							if reg in pams_CC:
								print change
								print 'shitting hell CC' 
								c+=1
								pam_not_destr +=1
								#exit()
						elif 'GG' in pam:
							reg = donor[donor.find(guide)+20:donor.find(guide)+23]
							if reg in pams_GG:
								print change
								print 'shitting hell GG'
								c+=1
								pam_not_destr +=1
								#exit()

						elif d > 6:
							print 'no good'
							c+= 1 
							in_but_far+=1
						else:
							good +=1
							g.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(change, position, pam, d, guide, donor))
					else:
						if d > 6:
							print d
							print 'not in guide but too far'
							not_in_but_far +=1
							c+= 1
							if d > 20:
								nibtf_20 +=1
						else:
							good+=1
							g.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(change, position, pam, d, guide, donor))
							
	print nb
	print 'total sequences {0} total filtered out {1} not_in_but_far (>6 <20) {2} not in donor but too far (>20) {3} in donor but too far {4} pam not destroyed {5} total files {6}'.format(total, c, not_in_but_far, nibtf_20, in_but_far, pam_not_destr, nb)
	print total - good - c


	for file in glob.glob('best_g+d_{0}_*_*_aa_filtered.txt'.format(nameee)):
		if os.stat(file).st_size > 0:
		   pass
		   #print "All good"
		else:
			os.remove(file)

	for file in glob.glob('best_g+d_{0}_*_*_aa.txt'.format(nameee)):
		os.remove(file)	

elif whattodo == 't':
	for file in glob.glob('best_g+d_*_*.txt'):
		if '_aa.txt' not in file:
			nb+= 1
			name = file.split('.')[0]
			nameee = file.split('_')[2]
			print name
			with open(file, 'r') as f:
				with open('{0}_filtered.txt'.format(name), 'w') as g:
					#g.write('seq/aa \tposition \tbest pam\t aminoacid dist from cut site \t\t guide \t\t\t donor\n')
					data = [i for i in f.read().split('\n') if len(i)>0]
					data = data[1:len(data)]
					#print len(data), data
		
					for line in data:
						total +=1 
						line = line.split('\t')
						change, position, pam, d, guide, donor = line[0], line[1], line[2], line[3], line[4], line[5]
						#print change, position, pam, d, guide, donor
						d = int(float(d))

						if guide in donor:
							print guide, donor 
							# fuck
							if 'CC' in pam:
								reg = donor[donor.find(guide)-3:donor.find(guide)]
								if reg in pams_CC:
									print change
									print 'shitting hell CC' 
									c+=1
									pam_not_destr +=1
									#exit()
							elif 'GG' in pam:
								reg = donor[donor.find(guide)+20:donor.find(guide)+23]
								if reg in pams_GG:
									print change
									print 'shitting hell GG'
									c+=1
									pam_not_destr +=1
									#exit()

							elif d > 6:
								print 'no good'
								c+= 1 
								in_but_far+=1
							else:
								good +=1
								g.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(change, position, pam, d, guide, donor))
						else:
							if d > 6:
								print d
								print 'not in guide but too far'
								not_in_but_far +=1
								c+= 1
								if d > 20:
									nibtf_20 +=1
							else:
								good+=1
								g.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(change, position, pam, d, guide, donor))
							
	print nb
	print 'total sequences {0} total filtered out {1} not_in_but_far (>6 <20) {2} not in donor but too far (>20) {3} in donor but too far {4} pam not destroyed {5} total files {6}'.format(total, c, not_in_but_far, nibtf_20, in_but_far, pam_not_destr, nb)
	print total - good - c


	for file in glob.glob('best_g+d_{0}_*_*_filtered.txt'.format(nameee)):
		if os.stat(file).st_size > 0 and '_aa.txt' not in file:
		   pass
		   #print "All good"
		else:
			os.remove(file)

	for file in glob.glob('best_g+d_{0}_*_*.txt'.format(nameee)):
		if '_aa.txt' not in file and '_filtered.txt' not in file:
			os.remove(file)		