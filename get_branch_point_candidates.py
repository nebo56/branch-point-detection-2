'''
Created on Jan 8, 2014

@author: Nejc Haberman

The script will filter out only sequences that ends with AG
'''


import sys

def filter_branch_point_seq (fin_sam, fout_sam):
    fin = open(fin_sam, "rt")
    fout = open(fout_sam, "w")
    line = fin.readline()
    while line[0] == '@':    #header of .SAM
        fout.write(line)
        line = fin.readline()
        
    while line:
        column = line.rstrip('\n').rstrip('\r').split('\t')
        if column[1] != '4':    #not error 
            seq = column[9]
            if column[1] == '0':    #same strand
                if str(seq[-2:]).upper() == "AG":
                    fout.write(line)
            if column[1] == '16':   #anti strand 
                if str(seq[0:2]).upper() == "CT":
                    fout.write(line)
        line = fin.readline()
    fin.close()
    fout.close()
    
if sys.argv.__len__() == 3:
    fin_sam = sys.argv[1]
    fout_sam = sys.argv[2]
    filter_branch_point_seq (fin_sam, fout_sam)
else:
    print "you need 2 arguments to run the script"
    quit()
