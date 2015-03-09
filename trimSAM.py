'''
Created on Jan 29, 2014

@author: Nejc Haberman


Script will trim reads that contains a genomic A mutation on first nucleotide  
'''

import sys


#method will return decreased column by 1. example: input XM:i:1, otput XM:i:0
def decrease_sam_column(col):
    tokens = col.rsplit(':')
    return tokens[0] + ':' + tokens[1] + ':' + str(int(tokens[-1]) - 1)


def trim_sam (fin_sam, fout_sam):
    finSam = open(fin_sam, "rt")
    foutSam = open(fout_sam, "w")
    line = finSam.readline()
    while line:
        if line[0] != '@':
            col = line.rstrip('\n').rsplit('\t')
            if col[1] != "4": #ignore unmapped reads
                seq = col[9]
                quality = col[10]
                missmatches = col[-2]
                tokens = missmatches.rsplit(':') #example: MD:Z:0A29
                if col[1] == "0": #same strand
                    if tokens[-1][0:2] == "0A": 
                        missmatches = missmatches.replace("0A","")
                        seq = seq[1:]
                        length = str(seq.__len__()) + "M"
                        quality = quality[1:]
                        XM = decrease_sam_column(col[-6])
                        NM = decrease_sam_column(col[-3])
                        position = str(int(col[3]) + 1)
                        foutSam.write(col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + position + '\t' + col[4] + '\t' + length + '\t' + col[6] + '\t' + col[7] + '\t' + col[8] + '\t' + seq + '\t' + quality + '\t' + col[11] + '\t' + col[12] + '\t' + XM + '\t' + col[14] + '\t' + col[15] + '\t' + NM + '\t' + missmatches + '\t' + col[18] + '\n')
                    else:
                        foutSam.write(line)
                        
                if col[1] == "16": #minus strand
                    if tokens[-1][-2:] == "T0": 
                        missmatches = missmatches.replace("T0","")
                        seq = seq[0:-1]
                        length = str(seq.__len__()) + "M"
                        quality = quality[0:-1]
                        XM = decrease_sam_column(col[-6])
                        NM = decrease_sam_column(col[-3])
                        position = str(int(col[3])) #position on the anti strand stays the same
                        foutSam.write(col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + position + '\t' + col[4] + '\t' + length + '\t' + col[6] + '\t' + col[7] + '\t' + col[8] + '\t' + seq + '\t' + quality + '\t' + col[11] + '\t' + col[12] + '\t' + XM + '\t' + col[14] + '\t' + col[15] + '\t' + NM + '\t' + missmatches + '\t' + col[18] + '\n')
                    else:
                        foutSam.write(line)
                        
        else:   #write header
            foutSam.write(line)
                
        line = finSam.readline()
    finSam.close()
    foutSam.close()

if sys.argv.__len__() == 3:
    fin_sam = sys.argv[1]
    fout_sam = sys.argv[2]
    trim_sam (fin_sam, fout_sam)
else:
    print "you need 2 arguments to run the script"
    quit()
