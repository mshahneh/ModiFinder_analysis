import pickle
from rdkit import Chem

import json
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import io
from IPython.display import SVG, display
import cairosvg
import math
import numpy as np
import uuid
import copy
import getpass
import sys
import os
import argparse
import requests
# add parent directory to path
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from SmallMol_Mod_Site_Localization import Compound_n as compound
from SmallMol_Mod_Site_Localization import ModificationSiteLocator as modSite
from SmallMol_Mod_Site_Localization import utils_n as utils
from SmallMol_Mod_Site_Localization import handle_network as handle_network
from SmallMol_Mod_Site_Localization import visualizer

def main(id_known, id_modified, method, args, output, data_path):
    try:
        with open(os.path.join(data_path, "cached_structures", id_known + ".pkl"), "rb") as f:
            mol_known = pickle.load(f)
        with open(os.path.join(data_path, "cached_compounds", id_known + ".pkl"), "rb") as f:
            value = pickle.load(f)
        known_compond = compound.Compound(value, mol_known, args={"should_fragment": False})
    except:
        known_compond = compound.Compound(id_known, args={"should_fragment": False})

    try:
        with open(os.path.join(data_path, "cached_structures", id_modified + ".pkl"), "rb") as f:
            mol_modified = pickle.load(f)
        with open(os.path.join(data_path, "cached_compounds", id_modified + ".pkl"), "rb") as f:
            value = pickle.load(f)
        modified_compond = compound.Compound(value, mol_modified, args={"should_fragment": False})
    except:
        modified_compond = compound.Compound(id_modified, args={"should_fragment": False})


    mod_site_obj = modSite.ModificationSiteLocator(known_compond, modified_compond)
    probabilities = mod_site_obj.generate_probabilities(method=method)
    try:
        probabilities = probabilities.tolist()
    except:
        try:
            for i in range(len(probabilities)):
                probabilities[i] = probabilities[i].tolist()
        except:
            pass
    graph = known_compond.distances.tolist()
    true_site = utils.calculateModificationSites(modified_compond.structure, known_compond.structure, False)[0]

    data = {"probabilities": probabilities, "graph": graph, "true_site": true_site}
    data["num_atoms"] = known_compond.structure.GetNumAtoms()
    # write data to json file
    with open(output, 'w') as outfile:
        json.dump(data, outfile)

    # url_bigger = "https://external.gnps2.org/gnpsspectrum?SpectrumID={}".format(args.id_modified)
    # value = requests.get(url_bigger).json()
    # mol_bigger = Chem.MolFromSmiles(value["annotations"][0]["Smiles"])

# if main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='run_random_experiment')
    parser.add_argument('--experiment_id', type=str, help='id of experiment to run')
    parser.add_argument('--experiments_metadata', type=str, help='path to csv file containing experiment metadata')
    parser.add_argument('--output', type=str, help='path to output file')
    parser.add_argument('--id_known', type=str, help='id of known compound')
    parser.add_argument('--id_modified', type=str, help='id of modified compound')
    args = parser.parse_args()

    experiments_metadata = pd.read_csv(args.experiments_metadata)
    # find experiment with id experiment_id
    experiments_metadata = experiments_metadata[experiments_metadata["id"] == args.experiment_id]
    # convert dataframe row 0 to dictionary
    experiments_metadata = experiments_metadata.to_dict(orient="records")[0]
    # add items of experiments_metadata to args
    args = vars(args)
    args.update(experiments_metadata)
    
    main(args["id_known"], args["id_modified"], experiments_metadata["method"], args, args["output"], "/home/user/LabData/Reza/data")

# python runner_random.py --experiment_id "2de2cb17df39a74931010eac88d56d7e4c63345c" --experiments_metadata "/home/user/Substructure_Assignment/my_implementation/SmallMol_Mod_Site_Localization_Analysis/results/experiments_meta_helpers.csv" --output "/home/user/Substructure_Assignment/my_implementation/SmallMol_Mod_Site_Localization_Analysis/results/test.json" --id_known "CCMSLIB00005464138" --id_modified "CCMSLIB00005464351"

