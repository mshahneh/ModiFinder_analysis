"""
This script is responsible for calculating the similar pairs (matched pairs of compounds) and their immidiate helpers.
"""
from urllib.request import urlopen
import pandas as pd
import os
import json
from tqdm import tqdm
import pickle
import re
import sys
# add parent directory to path
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from SmallMol_Mod_Site_Localization import utils_n as utils 

def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')


def get_elements(formula):
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    elements = [(element, count if count else 1) for element, count in elements]
    return elements

def get_diff(formula1, formula2):
    elements1 = get_elements(formula1)
    elements2 = get_elements(formula2)
    diff = {}
    for element, count in elements1:
        diff[element] = int(count)
    for element, count in elements2:
        if element in diff:
            diff[element] -= int(count)
        else:
            diff[element] = -int(count)
    
    # remove 0s
    diff = {key: value for key, value in diff.items() if value != 0}
    # dict to string
    diff = str(diff)
    return diff

def convert_to_formula(item):
    item = item.replace("{", "")
    item = item.replace("}", "")
    item = item.split(",")
    item = [re.sub(r'[\'\s]', '', i) for i in item]
    item = [i.split(":") for i in item]
    item = [[i[0], int(i[1])] for i in item]
    # item = sorted(item, key=lambda x: x[1], reverse=True)
    item = [i[0] + str(i[1]) for i in item]
    item = "".join(item)
    return item


def calculate_matches(library_compounds, library_structures, cirteria, difference_threshold_rate):
    """
    Calculates the matches between the library.
    It is assumed that the library is passed as a dictionary with unique key identifiers.
    Input:
        library_compounds: table with the compounds data
        library_structures: dictionary with the cached structures (built from SMILES strings using RDKit)
        difference_threshold_rate: the threshold for the difference in weight (in %)
    Output:
        table of matches
    """
    matches = pd.DataFrame()
    for i in tqdm(range(len(library_compounds)), ascii=True):
        compound1 = library_compounds.iloc[i]
        w1 = float(compound1['Precursor_MZ'])
        # get potential matches
        potential_matches = library_compounds[(library_compounds['Precursor_MZ'] > compound1['Precursor_MZ']) & (library_compounds['Precursor_MZ'] < compound1['Precursor_MZ'] * difference_threshold_rate)]
        for item in cirteria:
            potential_matches = potential_matches[potential_matches[item] == compound1[item]]
        for j in range(len(potential_matches)):
            compound2 = potential_matches.iloc[j]
            w2 = float(compound2['Precursor_MZ'])
            m1 = library_structures[compound1['spectrum_id']]
            m2 = library_structures[compound2['spectrum_id']]
            if m1.GetNumAtoms() < m2.GetNumAtoms() and m2.HasSubstructMatch(m1):
                numModificationSites = len(utils.calculateModificationSites(m2, m1))
                diff = convert_to_formula(get_diff(compound2['Formula_smiles'], compound1['Formula_smiles']))
                temp_data = {'id_bigger': compound2['spectrum_id'], 'id_smaller': compound1['spectrum_id'], 'num_modif_sites': numModificationSites, "diff_formula": diff}
                temp_data['weight_bigger'] = w2
                temp_data['weight_smaller'] = w1
                temp_data['difference'] = abs(w1 - w2)
                for item in cirteria:
                    item_lower = item.lower()
                    temp_data[item_lower] = compound1[item]
                matches = pd.concat([matches, pd.DataFrame(temp_data, index=[0])], ignore_index=True)
    
    return matches

def filter_helpers(id_target, weight_filter, pairs, eps = 0.5):
    potential_helpers_smaller = pairs[pairs['id_bigger'] == id_target]
    potential_helpers_smaller = potential_helpers_smaller[(potential_helpers_smaller['weight_smaller'] < weight_filter - eps) | (potential_helpers_smaller['weight_smaller'] > weight_filter + eps)]

    potential_helpers_bigger = pairs[pairs['id_smaller'] == id_target]
    potential_helpers_bigger = potential_helpers_bigger[(potential_helpers_bigger['weight_bigger'] < weight_filter - eps) | (potential_helpers_bigger['weight_bigger'] > weight_filter + eps)]

    id_helpers = list(potential_helpers_smaller['id_smaller'].values) + list(potential_helpers_bigger['id_bigger'].values)
    return id_helpers

def calculate_helpers(matches, output, max_num_modifications_allowed = 1):
    """
    Calculates the helpers for each compound.
    Input:
        matches: dataframe of matches
        max_num_modifications_allowed: maximum number of modification sites allowed
    Output:
        helpers: dictionary with the helpers, [key] = [list of helpers] where key is the compound id
    """
    potential_pairs = matches[(matches['num_modif_sites'] <= max_num_modifications_allowed)]
    helpers = {}
    for i, row in potential_pairs.iterrows():
        t_id_smaller = row["id_smaller"]
        t_weight_smaller = row["weight_smaller"]
        t_id_bigger = row["id_bigger"]
        t_weight_bigger = row["weight_bigger"]

        if helpers.get(t_id_smaller) is None:
            helpers[t_id_smaller] = []
        if helpers.get(t_id_bigger) is None:
            helpers[t_id_bigger] = []
        
        helpers_smaller = filter_helpers(t_id_smaller, t_weight_bigger, potential_pairs)
        if not os.path.exists(os.path.join(output, "helpers")):
            os.makedirs(os.path.join(output, "helpers"))
        with open(os.path.join(output, "helpers", "{}_{}.json".format(t_id_smaller, t_id_bigger)), 'w') as f:
            json.dump(helpers_smaller, f)
        
        helpers_bigger = filter_helpers(t_id_bigger, t_weight_smaller, potential_pairs)
        with open(os.path.join(output, "helpers", "{}_{}.json".format(t_id_bigger, t_id_smaller)), 'w') as f:
            json.dump(helpers_bigger, f)

def main(project_root, run_config):
    disable_rdkit_logging()

    data_folder = run_config["data_folder"]
    # if data_folder is not absolute path, make it absolute
    if not os.path.isabs(data_folder):
        data_folder = os.path.abspath(os.path.join(project_root, data_folder))

    libraries = pd.read_csv(os.path.join(project_root, "libraries.csv"))

    for i, row in libraries.iterrows():
        print(f"Processing library {row['Library']}")

        library_compounds_data = pd.read_csv(os.path.join(data_folder, "libraries", row['Library'] + "_compounds_data.csv"))

        library_compounds_data = library_compounds_data[library_compounds_data["Adduct"].isin(run_config['match_property']['supported_adducts'])]
        # add to meta how many pairs we loose
        cached_structures = pickle.load(open(os.path.join(data_folder, "cached_structures", row['Library'] + ".pkl"), "rb"))
        criteria = run_config['match_property']['matching_property_criteroa']
        difference_threshold_rate = run_config['match_property']['max_weight_diff_ratio']
        
        matches = calculate_matches(library_compounds_data, cached_structures, criteria, difference_threshold_rate)
        # add to meta how many pairs we loose
        if len(matches) == 0:
            print("No matches found for library: ", row['Library'])
            continue
        matches = matches[matches["num_modif_sites"] == 1]

        if not os.path.exists(os.path.join(data_folder, "matches")):
            os.makedirs(os.path.join(data_folder, "matches"))
        matches.to_csv(os.path.join(data_folder, "matches", row['Library'] + ".csv"), index=False)
        
        if len(matches) > 0:
            calculate_helpers(matches, data_folder, 1)

        print("Done", row['Library'], len(matches))

if __name__ == "__main__":
    project_root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    run_config = json.load(open(os.path.join(project_root, "run_config.json")))
    main(project_root, run_config)

