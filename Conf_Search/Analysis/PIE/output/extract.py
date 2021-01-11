import sys
import numpy as np
infile = sys.argv[1]
outfile = sys.argv[2]
output = open(outfile,'w')
f = open(infile,'r')
line = f.readline()
while line:
	line = line.split()
	if len(line)>2:
		if (line[0]=='ENERGY:'):
			#newline = line[1]+'\t'+line[6]+'\t'+line[7]+'\t'+str(float(line[6])+float(line[7]))+'\n'
			newline = line[11]+'\n'
			output.write(newline)
	line = f.readline()

output.close()
f.close()
