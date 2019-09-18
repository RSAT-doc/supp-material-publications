#!/usr/bin/env python
import sys
import re
file = sys.argv[1]
out = sys.argv[2]
print("Reading file: " + file)
with open (file, 'r') as varSeq, open (out + ".ref.fasta", 'w') as ref, open (out + ".alt.fasta", 'w') as alt:
	lines = [line for line in varSeq.readlines() if re.match('^\d', line)]
	for line in lines:
		cols = line.split('\t')
		header = ">" + cols[0] + ":" + cols[1]
		seq = cols[9].rstrip()
		if (cols[6] == cols[7]):
			#REFERENCIA
			print(header+"\n"+seq, file = ref)
			print("referencia")
		else:
			#ALTERNATIVO
			print(header+"\n"+seq, file = alt)
			print("alternativo")


	



	

