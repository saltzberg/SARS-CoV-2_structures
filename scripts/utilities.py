import Bio
import Bio.PDB
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from Bio.Blast import NCBIWWW


sequence_file = "../data/indexing/SARS_CoV_2.seq"
# Manually annotated PDB structures for SARS1 and SARS2
#
# Found via pblast of proteins in the PDB
#
### TWO ###

def create_sarscov_pdb_dictionaries(infile="../data/indexing/sars_cov_2_pdbblast.dat", s1={}, s2={}):
    f = open(infile,"r")

    seq_dict = parse_seq_file(sequence_file)

    for pn in get_protein_names("../data/indexing/protein_names.txt"):
        s1[pn] = []
        s2[pn] = []

    for l in f.readlines():
        fields = l.split(" ")
        pname = fields[0].split("_")[-1]
        pdbid = fields[1].split("_")[0]
        pdbchain = fields[1].split("_")[1]
        pct_id = float(fields[2])
        seqstart = int(fields[5]) 
        seqend = int(fields[6]) 
        offset = seqstart - int(fields[7]) 

        n_seq_res = len(seq_dict[pname])
        #if float(seqstart-seqend)/n_seq_res > 0.5 or seqstart-seqend > 20:
        #    continue            

        entry = (pdbid, pdbchain, seqstart, seqend, offset)

        # If pct_id > 98%, we assume this is a SARS-CoV-2 structure
        if pct_id > 98.0:
            s2[pname].append(entry)
        else:
            if pname != "S" and int(float(seqend-seqstart)/n_seq_res*100) > 10:
                s1[pname].append(entry)


    return s1, s2           


sars_cov_2_pdb_dictionary = {
    "nsp1" : [],
    "nsp2" : [],
    "nsp3" : [],
    "nsp4" : [],
    "nsp5" : [('5r7y',"a")],
    "nsp6" : [],
    "nsp7" : [('6m71', 'c'), ('7btf', 'c')],
    "nsp8" : [('6m71', 'b'), ('7btf', 'b')],
    "nsp9" : [('6w4b', 'a')],
    "nsp10" : [('6w4h', 'a')],
    "nsp11" : [],
    "nsp12" : [('6m71', 'a'), ('7btf', 'a')],
    "nsp13" : [],
    "nsp14" : [],
    "nsp15" : [('6vww', 'a')],
    "nsp16" : [('6w4h', 'a')],
    "E" : [],
    "M" : [],
    "N" : [('6m3m', 'a'), ('6vyo', 'a')],
    "orf3a" : [],
    "orf6" : [],
    "orf7a" : [],
    "orf8" : [],
    "orf9b" : [],
    "orf10" : [],
    "S" : [('6m17','e'),('6vsb','a'),('6vxx','a'),('6vyb','a')],
    "orf3b" : [],
    "orf7b" : [],
    "Protein14" : []
}

sars_cov_1_pdb_dictionary = {
    "nsp1" : [('2wct', "a"), ('3mj5', "a"), ('2gdt', "a")],
    "nsp2" : [],
    "nsp3" : [('5e6j', "a"), ('4mow', "a"), ('5tl6', "b"), ('3mj5', "a"), ('2fe8', "a"), ('5ye3', "a"), ('3e9s', "a")],
    "nsp4" : [],
    "nsp5" : [('6m03', 'a'), ('6nur', 'c'), ('1ysy', 'a')],
    "nsp6" : [],
    "nsp7" : [('2ahm','a')],
    "nsp8" : [('2ahm', 'b'), ('5f22', 'b'), ('6nur', 'b')],
    "nsp9" : [('1qz8', 'a'), ('1uw7', 'a'), ('3ee7', "a")],
    "nsp10" : [('5c8s', 'a'), ('2g9t', 'a'), ('5nfy', 'm'), ('3r24', 'b'), ('2fyg', 'a')],
    "nsp11" : [('3r23','b'), ('2g9t','a')],
    "nsp12" : [('6nur', 'a'), ('6nus','a')],
    "nsp13" : [('6jyt', 'a', 1)],
    "nsp14" : [('5c8s', 'b'),('5c8t', 'b'),  ('5nfy', "a")],
    "nsp15" : [('2h85', 'a'), ('2rhb', 'a'), ('2ozk', 'a')],
    "nsp16" : [('2xyr', "a"), ('2xyq', "a"), ('3r24', "a")],
    "E" : [('2mm4', 'a'), ('5x29', 'a')],
    "M" : [],
    "N" : [('2cjr', 'a'), ('2gib', 'a'), ('2ofz', 'a'), ('1ssk', 'a')],
    "orf3a" : [],
    "orf6" : [],
    "orf7a" : [('1xak', 'a'), ('1yo4', 'a')],
    "orf8" : [],
    "orf9b" : [('2cme', "a")],
    "orf10" : [],
    "S" : [('6acc', 'a'), ('5wrg', 'a'), ('5x58', 'a'), ('6nb6', 'a'), ('6crv', 'a')],
    "orf3b" : [],
    "orf7b" : [],
    "Protein14" : []
}

def get_protein_names():
    return sars_cov_2_pdb_dictionary.keys()

def get_first_residue_in_pdb(pdb, sequence):
    # Many homology models start from RESID=1, even if the first residues does not
    # correspond to residue = 1 in the sequence.  Find out what this is
    return 1

def is_subseq(x, y):
    it = iter(y)
    return all(c in it for c in x)

def get_sequence_from_pdb(pdb, nres="all"):
    parser = PDBParser(PERMISSIVE=1)
    struct = parser.get_structure("temp", pdb)
    seq = []
    resnums = []
    print(pdb)
    for r in struct.get_residues():
        #print(pdb, r.full_id[-1][1])
        resname = r.get_resname()
        rid = r.full_id[-1][1]
        seq.append(Bio.PDB.Polypeptide.three_to_one(resname))
        resnums.append(rid)
    first_res = min(resnums)
    last_res = max(resnums)
    fit_seq = ""

    for i in range(first_res, nres):
        if nres in resnums:
            fit_seq += seq[resnums.index(i)]
        else:
            break

    return fit_seq, first_res

def get_pdb_offset(pdb, sequence, nres=10):

    fit_seq, first_res = get_sequence_from_pdb(pdb, nres)

    if is_subseq(fit_seq, sequence):
        fits = []
        for i in range(len(sequence)):
            if sequence[i:i+len(fit_seq)] == fit_seq:
                fits.append(i)
        return fits

    else:
        raise("Exception: PDB sequence and target sequence do not match:", fit_seq)

def renumber_and_save_pdb(pdb, sequence, outpdb, diff=None):
    # Given a PDB file, return a 
    if diff is None:
        offsets = get_pdb_offset(pdb, sequence)
        if len(offsets) == 1:
            diff = offsets[0]
        else:
            print("Offset could be any of these:", offsets)
            print("Check manually and add to dictionary")
            exit()

    parser = PDBParser(PERMISSIVE=1)
    struct = parser.get_structure("temp", pdb)

    # Weird hack.  Need to change ID first to something ridiculous
    # then change it back. Otheriwse it yells at you that there are two
    # residues of the same ID.
    silly_large_integer = 89104
    # First chage by diff + a silly large integer
    for r in struct.get_residues():
        r.id = (r.full_id[3][0], r.full_id[3][1]+silly_large_integer+diff, r.full_id[3][2])
        new_tup = (r.full_id[0], r.full_id[1], r.full_id[2], r.id)
        r.full_id = new_tup
    # Then subtract silly large integer
    for r in struct.get_residues():
        r.id = (r.full_id[3][0], r.full_id[3][1]-silly_large_integer, r.full_id[3][2])
        new_tup = (r.full_id[0], r.full_id[1], r.full_id[2], r.id)
        r.full_id = new_tup
    
    io = PDBIO()
    io.set_structure(struct)
    io.save(outpdb)

def plot_coverage_single_protein(ax, sname, slen, s1_ranges, s2_ranges, 
                                label_offset=0):

    ax.set_frame_on(True)
    ax.set_ylim([0.0,0.9])
    ax.set_yticks([])
    ax.set_xlim([1,slen])
    ax.set_xticks([1,slen])
    ax.set_xticklabels([str(1),str(slen)], fontsize=14)
    #ax.text(slen/2, -0.2, s=sname, verticalalignment='bottom', horizontalalignment='center', fontstyle='italic', fontstretch='condensed', fontsize=14)
    ax.set_title(sname, fontsize=20)
    t = ax.text(-2, 0.15, s="SARS-CoV-2", verticalalignment='center', horizontalalignment='right', fontstyle='italic', fontstretch='condensed', fontsize=14)
    #t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))
    t = ax.text(-2, 0.45, s="SARS-CoV-1", verticalalignment='center', horizontalalignment='right', fontstyle='italic', fontstretch='condensed', fontsize=14)
    #t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))
    t = ax.text(-2, 0.75, s="Homology", verticalalignment='center', horizontalalignment='right', fontstyle='italic', fontstretch='condensed', fontsize=14)
    #t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))

    ax.plot((0,slen), (0.3,0.3), color=(0,0,0))
    ax.plot((0,slen), (0.6,0.6), color=(0,0,0))
    for rg in s2_ranges:
        ax.barh(0.15, 
                width=(rg[1]-rg[0]+1), 
                height=0.3, 
                left=(rg[0]-0.5+1), 
                color=(.16,.16,.88))

    for rg in s1_ranges:
        if rg[0] < 0:
            rg = (0, rg[1])
        if rg[1] < 0:
            rg = (0,1)
        ax.barh(0.45, 
                width=(rg[1]-rg[0]+1), 
                height=0.3, 
                left=(rg[0]-0.5+1) ,
                color=(0.4,0.4,0.8))

    return ax

def parse_seq_file(fname, protname=None):
    # Given a sequence file with lines:
    # >SARS_CoV_2_XXX
    # FASTASEQUENCE...
    # return a dictionary of protein names and sequences OR
    # return the sequence corresponding to protname
    seq_dict = {}
    f = open(fname, "r")
    sequence = ""
    for l in f.readlines():
        if len(seq_dict.keys()) == 0:
            pname = l.strip().split("_")[-1]
            seq_dict[pname] = sequence
        elif l[0] == ">":
            seq_dict[pname] = sequence
            sequence = ""
            pname = l.strip().split("_")[-1]
        else:
            sequence += l.strip()

    if protname is not None:
        return seq_dict[protname]
    else:
        return seq_dict

def get_protein_names(pname_file):
    f = open(pname_file, "r")
    protein_names = []
    f.readline() # Burn the first line
    for l in f.readlines():
        protein_names.append(l.strip())

    return protein_names

def get_ranges(ls):
    # for a list of numbers in order, return a list of ranges
    # i.e.
    #
    # [3,4,5,6,12,13,14,15,16]
    #
    # [(3,6),(12,16)]



    if len(ls) == 0:
        return []
    ls.sort()
    new_range_0 = ls[0]
    new_range_1 = ls[0]

    if len(ls) == 1:
        return [(new_range_0, new_range_1)]

    ranges = []

    for i in range(1, len(ls)):
        if ls[i] < new_range_1:
            raise Exception("list must be ordered", ls[i], "<", new_range_1)

        if i == len(ls)-1:
            ranges.append((new_range_0,ls[i]))
            break

        if ls[i] - new_range_1 == 1:
            new_range_1 = ls[i]
        else:
            ranges.append((new_range_0, new_range_1))
            new_range_0 = ls[i]
            new_range_1 = ls[i]

    return ranges
