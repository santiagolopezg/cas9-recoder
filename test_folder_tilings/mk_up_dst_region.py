# add flanking 100 nt up and downstream to sequence -> puts ATG at pos 100 using
# seq.find(ATG) and have no ATGs before that

import random
import math

n = 2
sites = ['ATG']
l = []

fasta_file = raw_input('Enter name of FASTA file: ')

seq = file.read(open(fasta_file, 'r')).split('\n')[1:]
seq = ''.join([i for i in seq if len(i) > 0])

p = seq.find('ATG')

stop_codons = ['TAA', 'TAG', 'TGA']

only_coding_seq = ''
c = ''
counter = p
while c not in stop_codons:
	
	c = seq[counter:counter+3]
	only_coding_seq += c
	#print c
	counter+=3


start = p
stop = counter

#print start, stop


upstream = p
downstream = int(math.fabs(counter - len(seq)))

# format upstream region
if upstream > 100:
	upstream_seq = seq[upstream-100:upstream]
elif upstream < 100 and upstream > 0:
	upstream_seq = (100-upstream)*'A'+seq[:upstream]
elif upstream == 0:
	upstream_seq = 100*'A'


if downstream == 0:
	downstream_seq = 100*'A'
elif downstream > 100:
	downstream_seq = seq[counter:counter+100]
elif downstream < 100 and downstream > 0:
	downstream_seq = seq[counter:]+(100-downstream)*'A'

#print upstream_seq, len(upstream_seq)
#print downstream_seq, len(downstream_seq)

ref_seq = upstream_seq + only_coding_seq + downstream_seq

name = fasta_file.split('.')[0]
with open(name+'.txt', 'w') as g:
	g.write(ref_seq)
	