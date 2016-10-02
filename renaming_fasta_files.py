#!/usr/bin/python
import os
from Bio import SeqIO



file_key = open('bifido_file_key.csv', "w")

count = 101


for filename in sorted(os.listdir(os.getcwd())):
    if filename.endswith(".fasta"):
        print filename
        
        file_key.write(filename + ',lac' + str(count) + '\n')
        
        renamed_fasta = open('lac%s.fasta' %(count), "w")

        orf_count = 1000
        
        for rec in SeqIO.parse(open(filename,"r"), "fasta"):
            if str(rec.seq).count("N") > 5:
                continue
            if len(str(rec.seq)) < 25: 
                continue

            renamed_fasta.write('>lac' + str(count) + '|' + str(count) + "_" + str(orf_count) + '\n'+str(rec.seq)+'\n')
            orf_count += 1


        count += 1

    else:
        continue


file_key.close()
    
