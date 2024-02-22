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
sys.path.insert(1, os.path.join(sys.path[0], ".."))

from SmallMol_Mod_Site_Localization import Compound_n as compound
from SmallMol_Mod_Site_Localization import ModificationSiteLocator as modSite
from SmallMol_Mod_Site_Localization import utils_n as utils
from SmallMol_Mod_Site_Localization import handle_network as handle_network
from SmallMol_Mod_Site_Localization import visualizer


def main(id_known, id_modified, method, args, output, data_path):
    # url_known = "https://external.gnps2.org/gnpsspectrum?SpectrumID={}".format(id_known)
    # value = requests.get(url_known).json()
    data = {}
    setting_args = {"ppm": args["mz_tol"], "fragmentation_depth": args["fragmentation_depth"]}
    try:
        with open(os.path.join(data_path, "cached_structures", id_known + ".pkl"), "rb") as f:
            mol_known = pickle.load(f)
        with open(os.path.join(data_path, "cached_compounds", id_known + ".pkl"), "rb") as f:
            value = pickle.load(f)
        known_compond = compound.Compound(value, mol_known, args=setting_args)
    except:
        known_compond = compound.Compound(id_known, args=setting_args)

    temp_args = copy.deepcopy(setting_args)
    if args["fragment_filter"] != "oracle" and args["fragment_filter"] != "oracle_nohelper":
        temp_args["should_fragment"] = False
    try:
        with open(os.path.join(data_path, "cached_structures", id_modified + ".pkl"), "rb") as f:
            mol_modified = pickle.load(f)
        with open(os.path.join(data_path, "cached_compounds", id_modified + ".pkl"), "rb") as f:
            value = pickle.load(f)
        modified_compond = compound.Compound(value, mol_modified, args=temp_args)
    except:
        modified_compond = compound.Compound(id_modified, args=temp_args)
        
    mod_site_obj = modSite.ModificationSiteLocator(known_compond, modified_compond, setting_args)
    true_site = utils.calculateModificationSites(
        modified_compond.structure, known_compond.structure, False
    )[0]
    shifted = [_[0] for _ in mod_site_obj.shifted]

    sirius_availability = False
    if args["fragment_filter"] == "msbuddy":
        mod_site_obj.main_compound.apply_msbuddy()
    if args["fragment_filter"] == "sirius" or args["fragment_filter"] == "oracle" or args["fragment_filter"] == "helpers" or args["fragment_filter"] == "oracle_nohelper":
        try:
            with open(os.path.join(data_path, "SIRIUS", id_known + "_fragmentationtree.json"), "rb") as f:
                sirius_data = json.load(f)
                mod_site_obj.main_compound.apply_sirius(sirius_data)
                sirius_availability = True
        except:
            print("no sirius data available")
            pass
    if args["fragment_filter"] == "helpers" or args["fragment_filter"] == "oracle":
        # try:
        
        with open(os.path.join(data_path, "helpers", "{}_{}.json".format(id_known, id_modified)), "rb") as f:
            helpers_data = json.load(f)

        helpers_array = []
        for helper in helpers_data:
            try:
                try:
                    with open(os.path.join(data_path, "cached_structures", helper + ".pkl"), "rb") as f:
                        mol_helper = pickle.load(f)
                    with open(os.path.join(data_path, "cached_compounds", helper + ".pkl"), "rb") as f:
                        value = pickle.load(f)
                    helper_compound = compound.Compound(value, mol_helper, args=setting_args)
                except:
                    helper_compound = compound.Compound(helper, args=setting_args)
                helpers_array.append(helper_compound)
            except:
                pass
        mod_site_obj.apply_helpers_compound_array(helpers_array, "intersection")
        data["helpers"] = helpers_data
        data["num_helpers"] = len(helpers_data)
        # except:
        #     data["helpers"] = []
        #     data["num_helpers"] = 0
        #     pass
    
    if args["fragment_filter"] == "oracle" or args["fragment_filter"] == "oracle_nohelper":
        mod_site_obj.main_compound.filter_fragments_by_atoms([true_site], shifted)

    shifted_only = True
    if args["use_unshifted"] == True or args["use_unshifted"] == "True":
        shifted_only = False
    
    probabilities = mod_site_obj.generate_probabilities(
        shifted_only=shifted_only, method=method
    ).tolist()
    graph = known_compond.distances.tolist()

    data["probabilities"] =  probabilities
    data["graph"] = graph
    data["true_site"] = true_site
    data["num_atoms"] = known_compond.structure.GetNumAtoms()
    data["cosine"] = mod_site_obj.cosine
    data["matched"] = len(mod_site_obj.matched_peaks)
    data["shifted"] = len(mod_site_obj.shifted)
    data["unshifted"] = len(mod_site_obj.unshifted)
    data["weight"] = known_compond.Precursor_MZ
    data["delta"] = abs(known_compond.Precursor_MZ - modified_compond.Precursor_MZ)

    data["peaks"] = len(known_compond.peaks)

    ambiguity, ratio = known_compond.calculate_peak_annotation_ambiguity(shifted)
    data["shifted_annotated_ratio"] = ratio
    data["shifted_annotated_ambiguity"] = ambiguity

    entropy = known_compond.calculate_annotation_entropy(shifted)
    data["shifted_annotated_entropy"] = entropy

    if args["fragment_filter"] == "sirius":
        data["sirius_availability"] = sirius_availability

    # write data to json file
    with open(output, "w") as outfile:
        json.dump(data, outfile)


    # url_bigger = "https://external.gnps2.org/gnpsspectrum?SpectrumID={}".format(args.id_modified)
    # value = requests.get(url_bigger).json()
    # mol_bigger = Chem.MolFromSmiles(value["annotations"][0]["Smiles"])


# if main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run_random_experiment")
    parser.add_argument("--experiment_id", type=str, help="id of experiment to run")
    parser.add_argument(
        "--experiments_metadata",
        type=str,
        help="path to csv file containing experiment metadata",
    )
    parser.add_argument("--output", type=str, help="path to output file")
    parser.add_argument("--id_known", type=str, help="id of known compound")
    parser.add_argument("--id_modified", type=str, help="id of modified compound")
    args = parser.parse_args()

    experiments_metadata = pd.read_csv(args.experiments_metadata)
    # find experiment with id experiment_id
    experiments_metadata = experiments_metadata[
        experiments_metadata["id"] == args.experiment_id
    ]
    # convert dataframe row 0 to dictionary
    experiments_metadata = experiments_metadata.to_dict(orient="records")[0]
    # add items of experiments_metadata to args
    args = vars(args)
    args.update(experiments_metadata)

    main(args["id_known"], args["id_modified"], experiments_metadata["method"], args, args["output"], "/home/user/LabData/Reza/data")