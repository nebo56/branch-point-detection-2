'''
Created on Aug 20, 2013

@author: Nejc

The script will flank the region in both directions
'''

import sys

def flank_positions(fin_fname, fout_fname, left_shift, right_shift):
    fin = open(fin_fname, "rt")
    fout = open(fout_fname, "w")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        pos1 = col[1]
        pos2 = col[2]
        pos1 = int(pos1) + int(left_shift)
        pos2 = int(pos2) + int(right_shift)
        if pos2 <= pos1:
		print "error: the intron is too short for flanking "
		print line
		pos1 = col[1]
		pos2 = col[2]
	fout.write(chr + '\t' + str(pos1) + '\t' + str(pos2) + '\n')
        line = fin.readline()

if sys.argv.__len__() == 5:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    left_shift = int(sys.argv[3])
    right_shift = int(sys.argv[4])
    flank_positions(fin_fname, fout_fname, left_shift, right_shift)
else:
    #print str(sys.argv.__len__())
    print "error:\t4 arguments are needed\n" + '\n' +"example:\t $ python bed_expand_positions.py input_fname.bed output_fname.bed left_shiftNUM right_shiftNUM"

