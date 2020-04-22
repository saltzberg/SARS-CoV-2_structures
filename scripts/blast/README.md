Scripts and database for finding SARS-CoV-2 and SARS-CoV-1 PDB structures:

1) Download the latest sequences from the PDB (~1min)

./ncbi-blast-2.10.0+/bin/update_blastdb.pl pdbaa; tar -zxvf pdbaa.tar.gz

2) Run blastp on the pdb database and save those >90% ID dat file

./ncbi-blast-2.10.0+/bin/blastp -db ./ncbi-blast-2.10.0+/pdbaa -query ../../data/indexing/SARS_CoV_2.seq -outfmt "6" | awk '$3>90.0{print $1, $2, $3, $4,"|", $7, $8, $9, $10}' > ../../data/indexing/sars_cov_2_pdbblast.dat

Output format:
Query_name PDB_Chain %ID align_length | seq_start seq_end pdb_start pdb_end
SARS_CoV_2_nsp1	2GDT_A	86.087	115	13	127	2	116
SARS_CoV_2_nsp1	1VS1_A	45.714	35	26	60	149	183

