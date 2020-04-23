# Create a figure showing the structure coverage for each chain in SARS-CoV-2

# Use the dictionaries in utilities.py
# sars_cov_2_pdb_dictionary - PDBID and chain IDs for SARS-CoV-2 structures
# sars_cov_1_pdb_dictionary - PDBID and chain IDs for SARS-CoV-1 structures

# Homology models not yet implemented.

import numpy
import glob
import Bio
import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
from matplotlib import pyplot as plt
from utilities import *

import warnings
warnings.filterwarnings("ignore")

# Directory setup
topdir = "../"
seq_file = topdir+"data/indexing/SARS_CoV_2.seq"

# First, get the all of the protein names
protein_names = get_protein_names(topdir+"data/indexing/protein_names.txt")

# Second, find the sequences associated with each protein. 
sd = parse_seq_file(seq_file)
protein_dict = {}
parser = PDBParser(PERMISSIVE=1)

# Make sars_cov_1/2 dictionaries
#pdbid_file = topdir+"data/indexing/sars_cov_2_pdbblast.dat" # used by default
sars_cov_1_pdb_dictionary, sars_cov_2_pdb_dictionary = create_sarscov_pdb_dictionaries()

# Loop over all proteins
for s in range(len(sd.keys())):

    fig, ax = plt.subplots(1, 1, figsize = (10,2))
    fig.subplots_adjust(left=-10)
    plt.tight_layout()
    sname = list(sd.keys())[s]

    seq_len = len(sd[sname])

    pdbs = sars_cov_2_pdb_dictionary[sname]

    sars2_resis = set()

    # Loop over all PDBs and get all residues covered
    for pdb in pdbs:
        pdb_file = topdir + "/pdbs/pdb"+pdb[0]+".ent"
        structure = parser.get_structure(pdb, pdb_file)
        chain_id = pdb[1]
        model = structure[0]
        resis = [r.get_id()[1]+pdb[-1] for r in model[chain_id].get_residues()]
        sars2_resis.update(resis)

    pdbs = sars_cov_1_pdb_dictionary[sname]

    sars1_resis = set()

    for pdb in pdbs:

        pdb_file = topdir + "/pdbs/pdb"+pdb[0]+".ent"
        structure = parser.get_structure(pdb, pdb_file)
        chain_id = pdb[1]
        model = structure[0]
        resis = [r.get_id()[1]+pdb[-1] for r in model[chain_id].get_residues()]

        sars1_resis.update(resis)



    # Plot coverage for SARS1 and SARS2
    # TODO - make this function universal.  Pass a single list of resis, y-value and name
    plot_coverage_single_protein(ax, sname, seq_len, get_ranges(list(sars1_resis)), get_ranges(list(sars2_resis)))

    plt.savefig("../figures/sequence_coverage/coverage_"+sname+".png", bbox_inches = "tight")

    plt.clf()
    plt.close(fig)
