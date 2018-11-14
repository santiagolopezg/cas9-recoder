#filter the filtered files for positions

import os, glob
import shutil

name = raw_input('Enter name of reference sequence file: ').split('.')[0]
whattodo = raw_input('Enter type of files to filter {r (recoding) // t (tiling)}: ')
positions = raw_input('Enter the positions of interest: ')

if whattodo == 'r':
	remove = raw_input('Do you wish to remove files with changes not ocurring in positions of interest? (yes / no) ')

if '-' in positions:
	poi = [i for i in xrange(int(positions.split('-')[0]),int(positions.split('-')[1])+1)]
else:
	poi = [int(i) for i in positions.split(',') if len(i)>0]

#poi = [706, 908, 909, 910, 911, 931, 932, 958, 959, 986, 990, 993, 703, 704]

files = []
#dir = 'g+d region of interest sgs1'

if whattodo == 'r':

	for file in glob.glob('best_g+d_*_*_*_aa_filtered.txt'):
		position = int(file.split('_')[4])
		aa = file.split('_')[3]
		#print position
		if position in poi:
			print position, aa
			files.append(file)
	print files
	
	if remove == 'yes':
		for f in glob.glob('best_g+d_{0}_*_*_aa_filtered.txt'.format(name)):
			if f not in files:
				os.remove(f)
	
			
if whattodo == 't':
	for file in glob.glob('best_g+d_{0}_*_filtered.txt'.format(name)):
		if '_aa_filtered.txt' not in file:
			nameme = file.split('.')[0]+'_for_pos.txt'
			with open(file, 'r') as f:
				good = []
				with open (nameme ,'w') as g:
					data = [i for i in f.read().split('\n') if len(i)>0]
					for line in data:
						pos = line.split('\t')[1]
						pos = pos.split(',')#.remove(',')
						pos = [int(e) for e in pos if len(e) > 0]
						#print pos, pos in poi
						for e in pos:
							print e, max(poi)
							if e in poi and max(poi)>=max(pos) and line not in good:
								print 'good'
								good.append(line)
					#exit()
					for i in good:
						print i
						g.write(i+'\n')
	



