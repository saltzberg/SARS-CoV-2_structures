import Bio
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from utilities import *
import glob
import os

# Using the SARS-COV-1 and SARS-CoV-2 dictionaries in utilities.py, 
# download each PDB file, extract the correct chain and place it in 
# the correct folder in ./proteins

class ChainSelect(Select):
    def __init__(self, chain):
        self.cid = chain
    def accept_chain(self, chain):
        if chain.id == self.cid:
            return True
        else:
            return False

def read_pdb_and_save_chain(pdbid, structure, chain, outdir):
    io = PDBIO()
    io.set_structure(structure)
    model = structure[0]
    ch = model[chain.upper()]
    outpdb_name = pdbid +"_"+chain+".pdb"
    io.save(outdir+outpdb_name, ChainSelect(chain.upper()))


topdir = "../"
seq_file = topdir+"indexing/SARS_CoV_2.seq"
protein_names = get_protein_names(topdir+"indexing/protein_names.txt")
sd = parse_seq_file(seq_file)

pdbl = PDBList()

parser = PDBParser(PERMISSIVE=1)

for n in protein_names:

    plist = sars_cov_2_pdb_dictionary[n]
    for p in range(len(plist)):
        pdb_id = plist[p][0]
        chain_id = plist[p][1]
        pdbl.retrieve_pdb_file(pdb_id, pdir=topdir+"pdbs/", file_format="pdb")
        pdb_file = topdir+"pdbs/pdb"+pdb_id+".ent"
        structure = parser.get_structure(pdb_id, pdb_file)
        outdir = topdir+"/proteins/"+n+"/SARS-CoV-2_pdb/"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        read_pdb_and_save_chain(pdb_id, structure, chain_id, outdir)

    plist = sars_cov_1_pdb_dictionary[n]
    for p in range(len(plist)):
        pdb_id = plist[p][0]
        chain_id = plist[p][1]
        pdbl.retrieve_pdb_file(pdb_id, pdir=topdir+"pdbs/", file_format="pdb")
        pdb_file = topdir+"pdbs/pdb"+pdb_id+".ent"
        structure = parser.get_structure(pdb_id, pdb_file)
        outdir = topdir+"/proteins/"+n+"/SARS-CoV-1_pdb/"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        read_pdb_and_save_chain(pdb_id, structure, chain_id, outdir)