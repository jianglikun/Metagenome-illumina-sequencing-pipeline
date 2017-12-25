#! /usr/bin/python

# -*- encoding: utf-8 -*-

import sys
import gzip
if __name__ == '__main__':
	mappedfile=sys.argv[1]
	rawfile=sys.argv[2]
	outfile=sys.argv[3]
	out_fh = open(outfile, 'w')
	mapped_reads = []
	with open(mappedfile,'r') as mp_fh:
		for line in mp_fh:
			mapped_reads.append(line.split("\t")[0])
	mapped_reads = set(mapped_reads)

	with  gzip.open(rawfile,'r') as rw_fh:
		flag = False
		for line in rw_fh:
			if line.startswith('@'):
				read=line.split(' ')[0].replace('@','')
				if read in mapped_reads:
					flag = False
					continue
				else:
					flag = True
			if flag:
				out_fh.write(line)
