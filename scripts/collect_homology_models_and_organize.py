# Script that collects homology models and organizes them into the appropriate folders
from bs4 import BeautifulSoup as BS 
import requests
import os
import zipfile
import gzip
import shutil
import glob
import sys
from utilities import *

databases = ["itasser", "Feiglab", "alphafold", "korkinlab", "swissprot", "modeller"]

# NCBI IDs for SARS-CoV-2 proteins

homology_dir = "../homology_datasets/"
model_dir = "../proteins/"
temp_dir = "./temp/"


def organize_models(pname_dict, prefix):
    for pname in pname_dict.keys():
        # Make sure we have a corresponding ncbi ID
        if pname_dict[pname] != "":
            prot_id = pname_dict[pname]
            pdb_file = homology_dir+prefix+"/"+prot_id+".pdb"

            # If we have a PDB file for this ID in the set, move
            # it to the correct homology model folder
            if os.path.exists(pdb_file):
                hmodel_dir = model_dir + pname + "/homology_models/"
                if not os.path.exists(hmodel_dir):
                    os.makedirs(hmodel_dir)
                shutil.copy(pdb_file, hmodel_dir + prefix +".pdb")

###############################
# I-Tasser
# "https://zhanglab.ccmb.med.umich.edu/C-I-TASSER/2019-nCov/"
#
# Zhang lab, U Michigan 
#
###############################

# For I-Tasser models
itasser_dict = {
    "nsp1" : "QHD43415_1",
    "nsp2" : "QHD43415_2",
    "nsp3" : "QHD43415_3",
    "nsp4" : "QHD43415_4",
    "nsp5" : "QHD43415_5",
    "nsp6" : "QHD43415_6",
    "nsp7" : "QHD43415_7",
    "nsp8" : "QHD43415_8",
    "nsp9" : "QHD43415_9",
    "nsp10" : "QHD43415_10",
    "nsp11" : "QHD43415_11",
    "nsp12" : "QHD43415_12",
    "nsp13" : "QHD43415_13",
    "nsp14" : "QHD43415_14",
    "nsp15" : "QHD43415_15",
    "orf3a" : "QHD43417",
    "E" : "QHD43418",
    "M" : "QHD43419",
    "orf6" : "QHD43420",
    "orf7a" : "QHD43421",
    "orf8" : "QHD43422",
    "N" : "QHD43423",
    "orf10" : "QHI42199",
    "S" : "QHD43416"
}

def update_ITasser_models(homology_dir=homology_dir, model_dir=model_dir, prefix="I-tasser"):
    # all models are combined in models.zip.
    # Each individual protein is named NCBI_ID.pdb.gz
    print("Downloading ITasser Models")
    url = "https://zhanglab.ccmb.med.umich.edu/C-I-TASSER/2019-nCov/models.zip"
    response = requests.get(url, stream=True)

    # Download models.zip
    with open(temp_dir+"models.zip",'wb') as f:
        for chunk in response.iter_content(chunk_size=1024): 
             # writing one chunk at a time to pdf file 
             if chunk: 
                 f.write(chunk) 

    # Unzip the main zipfile
    with zipfile.ZipFile(temp_dir+"models.zip", 'r') as zip_ref:
        zip_ref.extractall(temp_dir)

    # Now extract all of the individual PDBs
    for pgz in glob.glob(temp_dir+"models/*.pdb.gz"):
        with gzip.open(pgz, 'rb') as f_in:
            with open(pgz[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(pgz)

    if os.path.exists(homology_dir+prefix):
        shutil.rmtree(homology_dir+prefix)

    shutil.move(temp_dir+"models", homology_dir+prefix)

def organize_ITasser_models(homology_dir=homology_dir, model_dir=model_dir, prefix="I-tasser"):
    # Place homology models into the correct
    # Loop over all protein names
    print("Organizing ITasser Models")
    organize_models(itasser_dict, prefix)

################################
#
# SWISS-PROT
# "https://swissmodel.expasy.org/repository/species/2697049"
#
# 
#
################################

def update_swissprot_models(homology_dir=homology_dir, model_dir=model_dir, prefix="SWISS-PROT"):
    print("Downloading SWISSPROT models")
    url = "https://swissmodel.expasy.org/repository/species/2697049"

def organize_swissprot_models(homology_dir=homology_dir, model_dir=model_dir, prefix="SWISS-PROT"):
    print("Organizing SWISSPROT models")

################################
#
# AlphaFold
# 
#
# 
#
################################

def update_alphafold_models(homology_dir=homology_dir, model_dir=model_dir, prefix="SWISS-PROT"):
    print("Downloading AlphaFold models")
    url = ""

def organize_alphafold_models(homology_dir=homology_dir, model_dir=model_dir, prefix="SWISS-PROT"):
    print("Organizing AlphaFold models")
    pass

################################
#
# Korkin Lab
# 
#
# Worcester Polytech (WPI)
#
################################

def update_korkinlab_models(homology_dir=homology_dir, model_dir=model_dir, prefix="KorkinLab"):
    print("Downloading Korkin Lab models")
    #url = "https://swissmodel.expasy.org/repository/species/2697049"
    #response = requests.get(url, stream=True)

def organize_korkinlab_models(homology_dir=homology_dir, model_dir=model_dir, prefix="KorkinLab"):
    print("Organizing Korkin Lab models")
    pass

################################
#
# MODELLER
# Sali Lab
#
# UCSF
#
################################

def update_modeller_models(homology_dir=homology_dir, model_dir=model_dir, prefix="KorkinLab"):
    print("Downloading MODELLER models")
    pass

def organize_modeller_models(homology_dir=homology_dir, model_dir=model_dir, prefix="KorkinLab"):
    print("Organizing MODELLER models")
    pass

################################
#
# Feig Lab
# https://github.com/feiglab/sars-cov-2-proteins
#
# Worcester Polytech (WPI)
#
################################

feiglab_dict = {
    "nsp1" : "",
    "nsp2" : "nsp2",
    "nsp3" : "PL-PRO",
    "nsp4" : "nsp4",
    "nsp5" : "",
    "nsp6" : "nsp6",
    "nsp7" : "",
    "nsp8" : "",
    "nsp9" : "",
    "nsp10" : "",
    "nsp11" : "",
    "nsp12" : "",
    "nsp13" : "",
    "nsp14" : "",
    "nsp15" : "",
    "orf3a" : "ORF3a",
    "E" : "",
    "M" : "M_protein",
    "orf6" : "ORF6",
    "orf7a" : "",
    "orf7b" : "ORF7b",
    "orf8" : "ORF8",
    "N" : "",
    "orf10" : "ORF10",
    "S" : ""
}

def update_feiglab_models(homology_dir=homology_dir, model_dir=model_dir, prefix="FeigLab"):
    print("Downloading Feiglab models")
    url = "https://github.com/feiglab/sars-cov-2-proteins/archive/master.zip"
    response = requests.get(url, stream=True)
    with open(temp_dir+"master.zip",'wb') as f:
        for chunk in response.iter_content(chunk_size=1024): 
             # writing one chunk at a time to pdf file 
             if chunk: 
                 f.write(chunk) 

    # Unzip the main zipfile
    with zipfile.ZipFile(temp_dir+"master.zip", 'r') as zip_ref:
        zip_ref.extractall(temp_dir)

    if os.path.exists(homology_dir+prefix):
        shutil.rmtree(homology_dir+prefix)

    # Place the model files into the 
    shutil.move(temp_dir+"/sars-cov-2-proteins-master/FeigLab", homology_dir+prefix)

def organize_feiglab_models(homology_dir=homology_dir, model_dir=model_dir, prefix="FeigLab"):
    print("Organizing Feiglab models")
    organize_models(feiglab_dict, prefix)

###############
###############

def itasser(mode):
    if mode is "both":
        update_ITasser_models()
        organize_ITasser_models()
    elif mode is "organize":
        organize_ITasser_models()   
    else:
        update_ITasser_models() 

def swissprot(mode):
    if mode is "both":
        update_swissprot_models()
        organize_swissprot_models()
    elif mode is "organize":
        organize_swissprot_models()   
    else:
        update_swissprot_models() 

def korkinlab(mode):
    if mode is "both":
        update_korkinlab_models()
        organize_korkinlab_models()
    elif mode is "organize":
        organize_korkinlab_models()   
    else:
        update_korkinlab_models() 

def feiglab(mode):
    if mode is "both":
        update_feiglab_models()
        organize_feiglab_models()
    elif mode is "organize":
        organize_feiglab_models()   
    else:
        update_feiglab_models() 

def alphafold(mode):
    if mode is "both":
        update_alphafold_models()
        organize_alphafold_models()
    elif mode is "organize":
        organize_alphafold_models()   
    else:
        update_alphafold_models() 

def modeller(mode):
    if mode is "both":
        update_modeller_models()
        organize_modeller_models()
    elif mode is "organize":
        organize_modeller_models()   
    else:
        update_modeller_models() 

def run_all(mode):
    itasser(mode)
    swissprot(mode)
    korkinlab(mode)
    feiglab(mode)
    alphafold(mode)
    modeller(mode)

def print_help():
    print("Usage: collect_homology_models_and_organize.py mode database")
    print("--")
    print("mode")
    print("  h, help : display this help")
    print("  collect : (re)download all homology databases")
    print("  organize : organize currently downloaded homology databases")
    print("  both : download and organize homology databases")
    print("--")
    print("database")
    print("all : do for all databases")
    for db in databases:
        print(db + " :")

if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)

if len(sys.argv) > 1:
    mode = sys.argv[1]
else:
    mode = "both"

if len(sys.argv) > 2:
    database = sys.argv[2]
else:
    database = "all"


# Help
if mode == "h" or mode == "help":
    print_help()
    exit()
if mode=="both" or mode=="organize" or mode=="download":

    if database == "all":
        os.makedirs(temp_dir)
        run_all(mode)

    elif database in databases:
        # Make a temporary directory for downloads
        os.makedirs(temp_dir)
        if database =="itasser":
            itasser(mode)
        elif database == "modeller":
            modeller(mode)
        elif database=="feiglab":
            feiglab(mode)
        elif database=="swissprot":
            swissprot(mode)
        elif database=="korkinlab":
            korkinlab(mode)

        else:
            print("Database name is incorrect")
            print_help()
else:
    print("Mode is incorrect")
    print_help()




shutil.rmtree(temp_dir)








