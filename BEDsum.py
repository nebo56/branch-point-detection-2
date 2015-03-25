'''
Created on Mar 16, 2014

@author: Nejc Haberman


Script will sum together BED but BED needs to be sorted.
'''

import sys

def BEDsum(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()
    lastChr = None
    lastStart = None
    lastEnd = None
    lastStrand = None
    lastDistance = None
    cDNAsum = 0
    while line:
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        start = col[1]
        end = col[2]
        distance = col[4]
        cDNA = 1	# each read count as one
        strand = col[5]
        if lastChr == None or (lastChr == chr and lastStart == start and lastEnd == end and lastStrand == strand and lastDistance == distance):
            cDNAsum += int(cDNA)
        else:
            fout.write(lastChr + '\t' + lastStart + '\t' +  lastEnd + '\t' + str(cDNAsum) + '\t' + distance + '\t' + lastStrand + '\n')
            cDNAsum = int(cDNA)
        lastChr = chr
        lastStart = start
        lastEnd = end
        lastStrand = strand
        lastDistance = distance
        line = fin.readline()
    fout.write(lastChr + '\t' + lastStart + '\t' +  lastEnd + '\t' + str(cDNAsum) + '\t' + distance + '\t' + lastStrand + '\n')
    fout.close()
    fin.close()


if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    BEDsum(fname_in, fname_out)
else:
    print("python BEDsum.py <input_file> <output_file>")
