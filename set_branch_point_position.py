'''
Created on Jan 6, 2014

@author: skgthab
'''

'''
input: reads in bed format  with start and end positions
output: start position of the (branch point) read
'''

import sys

def set_branch_point_position (fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        if col[5] == '+':
            branch_point = int(col[1]) + 1
        else:
            branch_point = int(col[2]) - 1
        fout.write(col[0] + '\t' + str(branch_point) + '\t' + str(branch_point + 1) + '\t' + col[3] + '\t' + col[4] + '\t' + col[5] + '\n')
        line = fin.readline()
    fout.close()
    fin.close()


if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    set_branch_point_position (fname_in, fname_out)
else:
    print("python set_end_froset_branch_point_position.py <input_file> <output>")
