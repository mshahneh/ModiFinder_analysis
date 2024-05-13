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
import requests
from tqdm import tqdm
import copy
import getpass

import sys
import os
# add parent directory to path
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from SmallMol_Mod_Site_Localization import Compound_n as compound
from SmallMol_Mod_Site_Localization import ModificationSiteLocator as modSite
from SmallMol_Mod_Site_Localization import utils_n as utils
from SmallMol_Mod_Site_Localization import handle_network as handle_network
from SmallMol_Mod_Site_Localization import visualizer
from SmallMol_Mod_Site_Localization import calculate_scores_n as calculate_scores
import paper_figures.figures_handler as fh
from alignment import cosine_slow
from methods import *
from predict_modification_count import main as predict_modification_count
from rdkit import Chem

def disable_rdkit_logging():
        """
        Disables RDKit whiny logging.
        """
        import rdkit.rdBase as rkrb
        import rdkit.RDLogger as rkl

        logger = rkl.logger()
        logger.setLevel(rkl.ERROR)
        rkrb.DisableLog("rdApp.error")

disable_rdkit_logging()

dir_path = os.path.abspath('')
project_root = os.path.abspath(os.path.join(dir_path, os.pardir))

experiment_directory = os.path.join(project_root, 'experiments_settings', "experiments_meta_helpers.csv")
data_folder, results_directory, matches_directory, libraries, library_names = fh.get_basic_data(project_root)

libraries = list(library_names.keys())

res = {}

for library in libraries:
    with open(os.path.join(data_folder, "cached_structures", library + ".pkl"), "rb") as f:
        cached_structures = pickle.load(f)
    with open(os.path.join(data_folder, "cached_compounds", library + ".pkl"), "rb") as f:
        compounds_data = pickle.load(f)
    
    res[library] = {}

    pair_matches = pd.read_csv(library + ".csv")
    pair_matches = pair_matches[pair_matches['adduct'] == 'M+H']
    
    count = [0, 0, 0]

    false_positives = 0 # positive as in modiFinder works on it
    true_positives = 0
    true_negatives = 0
    false_negatives = 0

    count_bad = 0
    
    for index in tqdm(range(len(pair_matches))):
        try:
            id1 = pair_matches.iloc[index]["id_bigger"]
            id2 = pair_matches.iloc[index]["id_smaller"]
            args = {'ppm':40, 'mz_tolerance': 0.1}
            C1 = compound.Compound(compounds_data[id1], cached_structures[id1], args)
            C2 = compound.Compound(compounds_data[id2], cached_structures[id2], args)
            modiFinder = modSite.ModificationSiteLocator(C1, C2, args)

            if modiFinder.cosine <= 0.6:
                continue

            tuple1 = convert_to_SpectrumTuple(C1.peaks, C1.Precursor_MZ, C1.Charge)
            tuple2 = convert_to_SpectrumTuple(C2.peaks, C2.Precursor_MZ, C2.Charge)

            run_args = {
                'ppm': 40,
                'mz_tolerance': 0.1
            }
            pred = predict_modification_count(tuple1, tuple2, 'exhaustive', run_args=run_args)
            

            if pred:
                if pair_matches.iloc[index]["num_modif_sites"] == 1:
                    true_positives += 1
                else:
                    false_positives += 1
            else:
                if pair_matches.iloc[index]["num_modif_sites"] == 1:
                    false_negatives += 1
                else:
                    true_negatives += 1
        except Exception as e:
            count_bad += 1
            pass

    print("skipped: ", count_bad)
    res[library]["true_positives"] = true_positives
    res[library]["false_positives"] = false_positives
    res[library]["true_negatives"] = true_negatives
    res[library]["false_negatives"] = false_negatives
    res[library]['positives'] = true_positives + false_negatives
    res[library]['negatives'] = false_positives + true_negatives
    # print(count, count_bad, len(pair_matches[pair_matches["num_modif_sites"] <= 2]), len(pair_matches[pair_matches["num_modif_sites"] == 1]), len(pair_matches[pair_matches["num_modif_sites"] == 2]))
    # print(true_positives, false_positives, true_negatives, false_negatives)

    # print accuracy and precision and recall
    precision = true_positives / (true_positives + false_positives)
    recall = true_positives / (true_positives + false_negatives)
    accuracy = (true_positives + true_negatives) / (true_positives + true_negatives + false_positives + false_negatives)

    res[library]["precision"] = precision
    res[library]["recall"] = recall
    res[library]["accuracy"] = accuracy

    path = os.path.join(results_directory, "modification_count")
    if not os.path.exists(path):
        os.makedirs(path)
    with open(os.path.join(path, "results.json"), "w") as f:
        json.dump(res, f)
    print(library, precision, recall, accuracy)