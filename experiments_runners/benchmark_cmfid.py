# do the imports
import pickle
from rdkit import Chem
import argparse

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
# add parent directory to path
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from SmallMol_Mod_Site_Localization import Compound_n as compound
from SmallMol_Mod_Site_Localization import ModificationSiteLocator as modSite
from SmallMol_Mod_Site_Localization import utils_n as utils
from SmallMol_Mod_Site_Localization import handle_network as handle_network
from SmallMol_Mod_Site_Localization import visualizer

def create_link(library, id1, id2, smiles1, smiles2=None):
    base = "https://modsitelocalization.gnps2.org/"
    base = "http://localhost:5000/"
    usi1 = handle_network.generate_usi(id1, library)
    usi2 = handle_network.generate_usi(id2, library)
    if smiles2 is None:
        url = base + "?USI1=" + usi1 + "&USI2=" + usi2 + "&SMILES1=" + smiles1
    else:
        url = (
            base
            + "?USI1="
            + usi1
            + "&USI2="
            + usi2
            + "&SMILES1="
            + smiles1
            + "&SMILES2="
            + smiles2
        )
    return url

def task_compute_spectrum(lines):
    all_peaks = []
    current_energy = "NA"
    spectrum_identifier = 0
    # Reading Results
    for line in lines:
        if "energy" in line:
            current_energy = line.rstrip()
            spectrum_identifier += 1
            continue
        # if line is a new line, skip
        if line == "\n":
            break
        # Try to parse the peak intensities
        try:
            splits = line.rstrip().split(" ")
            mass = float(splits[0])
            intensity = float(splits[1])

            all_peaks.append((mass, intensity, current_energy, spectrum_identifier))
        except:
            pass

    mz = [peak[0] for peak in all_peaks]
    i = [peak[1] for peak in all_peaks]
    energy = [peak[2] for peak in all_peaks]
    identifier = [peak[3] for peak in all_peaks]

    peaks_df = pd.DataFrame()
    peaks_df['mz'] = mz
    peaks_df['i'] = i
    peaks_df['energy'] = energy
    peaks_df['identifier'] = identifier

    return peaks_df


def convert_to_conventional_peak(df):
    peaks = []
    # normalize the intensities
    max_i = max(df["i"])
    df["i"] = df["i"]/max_i
    for index, row in df.iterrows():
        peaks.append((row["mz"], row["i"]))
    return peaks

def peak_score_value(peak1, peak2, Precursor_MZ, Charge, method="cosine", use_intensity=True, fragment_mz_tolerance=0.1, fragment_ppm_tolerance=10):
    if not use_intensity:
        # copy the peaks
        peak1 = peak1.copy()
        peak2 = peak2.copy()
        for i in range(len(peak1)):
            peak1[i] = (peak1[i][0], 1)
        for i in range(len(peak2)):
            peak2[i] = (peak2[i][0], 1)
    
    if method == "cosine":
        from SmallMol_Mod_Site_Localization.alignment_n import _cosine_fast
        main_spectrum = utils.convert_to_SpectrumTuple(peak1, Precursor_MZ, Charge)
    modified_spectrum = utils.convert_to_SpectrumTuple(peak2, Precursor_MZ, Charge)
    cosine, matched_peaks = _cosine_fast(
        main_spectrum, modified_spectrum, fragment_mz_tolerance, fragment_ppm_tolerance, False
    )
    return cosine, matched_peaks

def calculate_probabilities(main_compound, mod_compound, filename, data_folder, use_intensity=True, fragment_mz_tolerance=0.1, fragment_ppm_tolerance=40):
    smiles = []
    with open(os.path.join(data_folder, "cfmid_exp/synth_output",filename), "r") as f:
        cfmid_keys = f.readlines()
        for cfmid_key in cfmid_keys:
            # remove the \n
            cfmid_key = cfmid_key[:-1]
            if len(cfmid_key) > 3 and len(cfmid_key.split(" ")) == 2:
                smiles.append(cfmid_key.split(" ")[1])
    
    print(os.path.join(data_folder, "cfmid_exp/synth_output",filename))
    print(smiles)
    
    preds = []
    with open(os.path.join(data_folder, "cfmid_exp/cfmid_preds", filename), "r") as f:
        all_predictions = f.readlines()
        # section the file into chunks when #In-silico ESI-MS/MS [M+H]+ Spectra is found
        # each chunk is a different prediction
        chunk = []
        for line in all_predictions:
            if line.startswith("#In-silico"):
                if len(chunk) > 0:
                    preds.append(chunk)
                chunk = [line]
            else:
                chunk.append(line)
        preds.append(chunk)


    for i in range(len(preds)):
        preds[i] = task_compute_spectrum(preds[i])

    probabilities = [0 for i in range(main_compound.structure.GetNumAtoms())]
    for i in range(len(preds)):
        predicted_modified_mol = Chem.MolFromSmiles(smiles[i])
        modif_site = utils.calculateModificationSites(predicted_modified_mol, main_compound.structure, False)[0]
        cosine, matched_peaks = peak_score_value(main_compound.peaks, convert_to_conventional_peak(preds[i]), main_compound.Precursor_MZ, main_compound.Charge, use_intensity=use_intensity, fragment_mz_tolerance=fragment_mz_tolerance, fragment_ppm_tolerance=fragment_ppm_tolerance)
        probabilities[modif_site] = cosine

    return probabilities

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_path", type=str)
    args = parser.parse_args()
    data_folder = args.data_path
    columns = [
        "id_smaller",
        "id_bigger",
    ]
    scoring_systems = ["is_max", "dist_from_max", "average_dist_from_max", "average_dist", "regulated_exp", "ranking_loss"]
    new_columns = []
    for scoring_system in scoring_systems:
        new_columns.append(scoring_system)
    columns += new_columns

    df = pd.DataFrame(columns=columns)
    args = {"filter_peaks_method": "intensity", "filter_peaks_variable": 0.01, "mz_tolerance": 0.05, "ppm":40}
    file_names = os.listdir(os.path.join(data_folder, "cfmid_exp/cfmid_preds"))

    for filename in file_names:
        try:
            id_smaller, id_bigger = filename.split(".")[0].split("_")
            with open(os.path.join(data_folder, "cached_compounds", id_smaller + ".pkl"), "rb") as f:
                smaller_compound = pickle.load(f)
            with open(os.path.join(data_folder, "cached_structures", id_smaller + ".pkl"), "rb") as f:
                smaller_struct = pickle.load(f)
            with open(os.path.join(data_folder, "cached_compounds", id_bigger + ".pkl"), "rb") as f:
                bigger_compound = pickle.load(f)
            with open(os.path.join(data_folder, "cached_structures", id_bigger + ".pkl"), "rb") as f:
                bigger_struct = pickle.load(f)
            main_compound = compound.Compound(smaller_compound, smaller_struct, args)
            mod_compound  = compound.Compound(bigger_compound, bigger_struct, args)
            modSiteLocator = modSite.ModificationSiteLocator(main_compound, mod_compound, args)
            true_modification_site = utils.calculateModificationSites(mod_compound.structure, main_compound.structure, False)[0]

            probabilities = calculate_probabilities(main_compound, mod_compound, filename, data_folder)

            val = {}
            val["id_smaller"] = id_smaller
            val["id_bigger"] = id_bigger

            for scoring_system in scoring_systems:
                val[scoring_system] = modSiteLocator.calculate_score(true_modification_site, scoring_system, probabilities)

            for col in columns:
                    if col not in val:
                        val[col] = None
            df = pd.concat(
                [
                    df,
                    pd.DataFrame(
                        val,
                        index=[0],
                    ),
                ]
            )

            df.to_csv("cfmid_scores.csv")
        except:
            pass
