#!/bin/bash

# Scripts and database for finding SARS-CoV-2 and SARS-CoV-1 PDB structures:
# 1) Download the latest sequences from the PDB (~1min)

./update_blastdb.pl pdbaa; tar -zxvf pdbaa.tar.gz

#2) Run blastp on the pdb database and save those >90% ID dat file (15s)

./blastp -db ./pdbaa -query ../../data/indexing/SARS_CoV_2.seq -outfmt "6" | awk '$3>90.0{print $1, $2, $3, $4,"|", $7, $8, $9, $10}' > ../../data/indexing/sars_cov_2_pdbblast.dat

