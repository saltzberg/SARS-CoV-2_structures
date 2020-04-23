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

def read_pdb_and_save_chain_with_offset(pdbid, structure, chain, outdir, offset=0):
    io = PDBIO()
    io.set_structure(structure)
    model = structure[0]
    ch = model[chain.upper()]
    outpdb_name = pdbid +"_"+chain+".pdb"

    silly_large_integer = 89104
    # First change by offset + a silly large integer
    for r in structure.get_residues():
        r.id = (r.full_id[3][0], r.full_id[3][1]+silly_large_integer+offset, r.full_id[3][2])
        new_tup = (r.full_id[0], r.full_id[1], r.full_id[2], r.id)
        r.full_id = new_tup
    # Then subtract silly large integer
    for r in structure.get_residues():
        r.id = (r.full_id[3][0], r.full_id[3][1]-silly_large_integer, r.full_id[3][2])
        new_tup = (r.full_id[0], r.full_id[1], r.full_id[2], r.id)
        r.full_id = new_tup

    io.save(outdir+outpdb_name, ChainSelect(chain.upper()))


topdir = "../"
seq_file = topdir+"data/indexing/SARS_CoV_2.seq"


# Make sars_cov_1/2 dictionaries
#pdbid_file = topdir+"data/indexing/sars_cov_2_pdbblast.dat" # used by default
sars_cov_1_pdb_dictionary, sars_cov_2_pdb_dictionary = create_sarscov_pdb_dictionaries()

protein_names = get_protein_names(topdir+"data/indexing/protein_names.txt")
sd = parse_seq_file(seq_file)
pdbl = PDBList()

parser = PDBParser(PERMISSIVE=1)

# Clear all existing pdbs
os.system("rm -rf topdir/proteins/*/SARS*/*")

for n in protein_names:

    plist = sars_cov_2_pdb_dictionary[n]
    for p in range(len(plist)):
        pdb_id = plist[p][0]
        chain_id = plist[p][1]
        offset = plist[p][-1]
        pdbl.retrieve_pdb_file(pdb_id, pdir=topdir+"pdbs/", file_format="pdb")
        pdb_file = topdir+"pdbs/pdb"+pdb_id+".ent"
        structure = parser.get_structure(pdb_id, pdb_file)
        outdir = topdir+"/proteins/"+n+"/SARS-CoV-2_pdb/"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        read_pdb_and_save_chain_with_offset(pdb_id, structure, chain_id, outdir, offset)

    plist = sars_cov_1_pdb_dictionary[n]
    for p in range(len(plist)):
        pdb_id = plist[p][0]
        chain_id = plist[p][1]
        offset = plist[p][-1]
        pdbl.retrieve_pdb_file(pdb_id, pdir=topdir+"pdbs/", file_format="pdb")
        pdb_file = topdir+"pdbs/pdb"+pdb_id+".ent"
        structure = parser.get_structure(pdb_id, pdb_file)
        outdir = topdir+"/proteins/"+n+"/SARS-CoV-1_pdb/"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        read_pdb_and_save_chain_with_offset(pdb_id, structure, chain_id, outdir, offset)


# TODO:  Add re-refined structures of Tristan Croll and ThornLab
thornlab = "https://github.com/thorn-lab/coronavirus_structural_task_force"