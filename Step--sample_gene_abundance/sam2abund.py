#! /usr/bin/env python

import getopt, sys
import redis

from collections import Counter

def ras_line(l):
  return l.strip().split("\t")

def CIGAR_parser(s):
    d = {}
    i = 0
    while i < len(s):
        for j in range(i+1, len(s)):
            if s[j] in ('M','I','D','N','S','H','P','=','X'):
                if s[j] in d:
                    d[s[j]] +=int(s[i:j])
                else:
                    d[s[j]] = int(s[i:j])
                i = j+1
                break

    return d


def sam2count(samfile, minid):
    rc = Counter()
    with open(samfile, "r") as fi:
        for line in fi:
            o = ras_line(line)
            if o[0][0] !="@" and o[2][-1] != "*" and o[5][-1] != "*":
                cigar = CIGAR_parser(o[5])
                if cigar['M'] >=int(minid*sum(cigar.values())):
                    rc[o[2]] +=0.5
    rcf = {}
    for gene in rc:
        if rc[gene]>=1: #filter out low confident hits
            rcf[gene] = rc[gene]
    #print(gene + '\t' + str(rcf[gene]))
    return rcf


def count2rpk(rc, lenprof):
    rpk = {}
    for gene in rc.keys():
        v = rc[gene]/float(lenprof[gene])
        rpk[gene]=v*1000 #per kilobase
    return rpk

def rpk2tpm(rpk):
    totalRPK =sum(rpk.values())

    totalRPK = totalRPK/1000000.0# normalize to 1M to get the scaling factor

    abunTable = {}
    for gene in rpk:
        abunTable[gene] = rpk[gene]/totalRPK

    
    return abunTable

def readGeneLength(genelengthfile):
    res = {}
    with open(genelengthfile,'r') as fh:
        for line in fh:
            line = line.strip()
            items = line.split('\t')
            res[items[0]] = items[1]
    return res

def main():
    #-----parse arguments from inputs------------
    try:
        options, remainder = getopt.getopt(sys.argv[1:],"",  ['sam=','genelength=', 'min_identity='])
    except getopt.GetoptError as err:
        print (err)
        sys.exit(2)

    samfile=''
    minIdentity = 0.95

    for op,value in options:
            if op=='--sam':
                    samfile=value
            elif op=='--genelength':
                genelengthfile = value
            elif op=='--min_identity':
                minIdentity = float(value)
            else:
                    print("sam2qha transforms a sam file produced by aligners to a binary abundance file.\nUsage: python sam2abund.py --sam samfile  --genelength  --min_identity required minimum identity (0-1) of a read to the reference, default 0.95")
                    sys.exit()



    #---------------sam 2 readcount------------
    readCount = sam2count(samfile, minIdentity)
    

    print("Successfully calculated read counts.")

    #---------fetch gene length from Redis server------------

    geneLengths = readGeneLength(genelengthfile)
    print("Successfully fetched gene length profile from redis.")

    #--------------readcount 2 abundance--------------
    RPK = count2rpk(readCount, geneLengths)
    abunTable = rpk2tpm(RPK)
    print("Successfully calculated relative abundance for each gene.")
    fh1 = open(samfile.split('.')[0] + '.readCount.txt','w')
    fh2 = open(samfile.split('.')[0] + '.tpm.txt','w')
    genes = readCount.keys()
    for gene in genes:
        fh1.write(gene + '\t' + str(readCount[gene]) + '\n')
        fh2.write(gene + '\t' + str(100*abunTable[gene]) + '\n')
    fh1.close()
    fh2.close()


if __name__=="__main__":
	main()


#------------------------
