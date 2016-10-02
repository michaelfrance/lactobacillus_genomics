#!/usr/bin/python
from Bio import SeqIO                                                               
import sys                                                                          
import fileinput
import shutil
import os
import sys
from Bio import AlignIO
from Bio.Alphabet import generic_protein
from Bio.Nexus import Nexus
from Bio.Align.Applications import ClustalwCommandline
import time

paralog = 0
core = 0
accessory = 0
single = 0

cog_fixed = open('crisp_access_for_cog.fasta','w')
single_cog_fixed = open('crisp_singlecog.fasta',"w")

crisp_ids = ["lac102","lac103","lac104","lac105","lac106","lac107","lac108","lac109","lac110","lac111","lac112","lac113","lac114","lac115","lac116"]

iners_ids = ["lac120","lac121","lac122","lac123","lac124","lac125","lac126","lac127","lac128","lac129","lac130","lac131","lac132","lac133","lac134"]


for line in fileinput.input(sys.argv[1]):
    line = line.strip()
    orf = line.split(" ")[0]
    orf = orf.replace(':','')
    #print orf

    protein_ids_init = line.split(" ")
    protein_ids_init.pop(0)
    protein_ids = list()
    
    for protein_id in protein_ids_init:
        protein_id_2 = protein_id.split("|")[0]

        if protein_id_2 in crisp_ids:
            protein_ids.append(protein_id)
    #print protein_ids
    #time.sleep(5)


    count=0
    protein_org_set = list()
    for i in protein_ids:
       protein_org = i.split("|")[1]       
       protein_org = protein_org.split("_")[0]
       protein_org_set.append(protein_org)
       #print protein_org_set 
       count +=1


    size1=len(protein_org_set)
    protein_org_set = list(set(protein_org_set))
    #print protein_org_set
    size2=len(protein_org_set)
 
    #if size1 > size2:
    #    paralog += 0    
    if size2 == 15 or size2 == 14:
        core += 1
    elif 1 < size2 < 14:
        accessory += 1

        orf_seq = open('%s.fasta' % (orf),'w')
        seqiter = SeqIO.parse(open(sys.argv[2]), 'fasta')                                    
        SeqIO.write((seq for seq in seqiter if seq.id in protein_ids), orf_seq, "fasta")
        orf_seq.close()
        filename = orf + '.fasta'

        #this loops opens the resulting multifasta from the previous step and formats the names a bit again and then writes the information to a multifasta files with data from each cluster
        count = 0
        for rec in SeqIO.parse(open(filename,"r"), "fasta"):            
            if count == 0:
                rec.id = rec.id.split("|")[1]
                rec.id = rec.id.split("_")[0]
                rec.id = rec.id.replace("-","")

                cog_fixed.write('>' + orf + '_' + rec.id +'\n'+str(rec.seq)+'\n')
                count += 1

            
    elif size2 == 1:
        single += 1
        orf_seq = open('%s.fasta' % (orf),'w')
        seqiter = SeqIO.parse(open(sys.argv[2]), 'fasta')                                    
        SeqIO.write((seq for seq in seqiter if seq.id in protein_ids), orf_seq, "fasta")
        orf_seq.close()
        filename = orf + '.fasta'

        #this loops opens the resulting multifasta from the previous step and formats the names a bit again and then writes the information to a multifasta files with data from each cluster
        count = 0
        for rec in SeqIO.parse(open(filename,"r"), "fasta"):            
            if count == 0:
                rec.id = rec.id.split("|")[1]
                rec.id = rec.id.split("_")[0]
                rec.id = rec.id.replace("-","")

                single_cog_fixed.write('>' + orf + '_' + rec.id +'\n'+str(rec.seq)+'\n')
                count += 1





cog_fixed.close()
single_cog_fixed.close()

print("paralog\t" + str(paralog))
print("core\t" + str(core))
print("accessory\t" + str(accessory))
print("single\t" + str(single))


