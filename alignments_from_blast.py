#!/usr/bin/python
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio.Align.Applications import ClustalwCommandline
import sys



gene_name = sys.argv[1].split(".")[0]

fasta_out = open('%s.fasta' %(gene_name), "w")

fasta_name = gene_name + '.fasta'

blast_records = NCBIXML.parse(open(sys.argv[1]))
count = 0
for blast_record in blast_records:
        for alignment in blast_record.alignments:
            sequence_name = alignment.hit_def
            count += 1
            if count < 26:
                counter2 = 0
                for hsp in alignment.hsps:
                     if counter2 < 1:

                     	 #sequence_name = alignment.title.split("|")[2]
            	 	
               	         fasta_out.write('>' + str(sequence_name) + '\n' + hsp.sbjct + '\n')
                         counter2 += 1

cline = ClustalwCommandline("clustalw2", infile=fasta_name, output="PHYLIP")
cline()
