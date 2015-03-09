'''
Created on Mar 5, 2014

@author: Nejc Haberman

Script will count a genomic transitions on first nucleotide from each read. Input file must be in SAM format.

'''

import sys

def count_transitions (fin_sam, fout_log):
    fin = open(fin_sam, "rt")
    fout = open(fout_log, "w")
    sam = fin.readline()
    reads = {}
    A,C,T,G = 0,0,0,0
    Ac,Cc,Tc,Gc = 0,0,0,0 #sum of all acids on 1. nt 
    read_count = 0
    while sam[0] == '@':    #header of .SAM
        sam = fin.readline()
    while sam:
        tokens = sam.rstrip('\n').rstrip('\r').split('\t')
        flag = tokens[1]
        if flag != "4":
            seq = tokens[9]
            transitions = tokens[-2]    #transition column is always one before the last column

            #counting genomic transitions
            if flag == "0": #same strand
                A += transitions.count(':0A')
                C += transitions.count(':0C')
                T += transitions.count(':0T')
                G += transitions.count(':0G')
                 
                #counting all transitions
                Ac += transitions.count('A')
                Cc += transitions.count('C')
                Tc += transitions.count('T')
                Gc += transitions.count('G')
                
            if flag == "16": #anti strand 
                A += transitions.count(':0T')
                C += transitions.count(':0G')
                T += transitions.count(':0A')
                G += transitions.count(':0C')
                
                #counting all transitions
                Ac += transitions.count('A')
                Cc += transitions.count('C')
                Tc += transitions.count('T')
                Gc += transitions.count('G')

            read_count += 1
        sam = fin.readline()

    # write results
    sum = A + C + T + G #sum of reads with first transition
    pA = (float(A)/float(sum)) * 100.0
    pC = (float(C)/float(sum)) * 100.0
    pT = (float(T)/float(sum)) * 100.0
    pG = (float(G)/float(sum)) * 100.0
    ratio_to_total = (float(sum)/float(read_count)) * 100.0
    fout.write("Number of first nt transitions on genomic: \n")
    fout.write("A: " + str(A) + "\t(" + str(round(pA,2)) + " %)\n")
    fout.write("C: " + str(C) + "\t(" + str(round(pC,2)) + " %)\n")
    fout.write("T: " + str(T) + "\t(" + str(round(pT,2)) + " %)\n")
    fout.write("G: " + str(G) + "\t(" + str(round(pG,2)) + " %)\n")
    fout.write("number of reads with first transition: " + str(sum) + "\t(" + str(round(ratio_to_total,2)) + " %)\n")

    sumC = Ac + Cc + Tc + Gc #sum of rall transitions
    pAc = (float(Ac)/float(sumC)) * 100.0
    pCc = (float(Cc)/float(sumC)) * 100.0
    pTc = (float(Tc)/float(sumC)) * 100.0
    pGc = (float(Gc)/float(sumC)) * 100.0
    fout.write("\nNumber of all transitions on genomic: \n")
    fout.write("A: " + str(Ac) + " \t(" + str(round(pAc,2)) + " %)\n")
    fout.write("C: " + str(Cc) + " \t(" + str(round(pCc,2)) + " %)\n")
    fout.write("T: " + str(Tc) + " \t(" + str(round(pTc,2)) + " %)\n")
    fout.write("G: " + str(Gc) + " \t(" + str(round(pGc,2)) + " %)\n")
    ratio_to_totalC = (float(sumC)/float(read_count)) * 100.0
    fout.write("number of all transitions: " + str(sumC) + "\t(" + str(round(ratio_to_totalC,2)) + " %)\n")

if sys.argv.__len__() == 3:
    fin_sam = sys.argv[1]
    fout_log = sys.argv[2]
    count_transitions(fin_sam, fout_log)
else:
    print "you need 2 arguments to run the script"
    quit()


        
