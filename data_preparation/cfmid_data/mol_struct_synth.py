"""
This code creates the permutations of the modification around the molecule.
"""

import sys
import os
import argparse
from rdkit import Chem
import pickle
import pandas as pd
import requests
sys.path.insert(1, os.path.join(sys.path[0], "../.."))
from SmallMol_Mod_Site_Localization import utils_n as utils
import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl
import math

def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog("rdApp.error")

def get_data(matches_path, batch_number, number_of_batches):
    """
    Get the data for the current batch
    """
    matches = pd.read_csv(matches_path)
    batch_size = math.floor(len(matches) / number_of_batches)
    if batch_size == 0:
        batch_size = 1
    if batch_number == number_of_batches - 1:
        batch = matches.iloc[batch_number * batch_size :]
    else:
        batch = matches.iloc[batch_number * batch_size : (batch_number + 1) * batch_size]
    return batch

def create_permutations(id_smaller, id_bigger, cached_structures, output_dir):
    """
    Create the structures
    """
    mol_smaller = cached_structures[id_smaller]
    mol_bigger = cached_structures[id_bigger]
    structs = utils.generate_possible_stuctures(mol_bigger, mol_smaller)

    df = pd.DataFrame(columns=["atom_index", "struct_smiles"])
    for i, struct in enumerate(structs):
        atom_index = struct[0]
        mol_struct = struct[1]
        df.loc[i] = [atom_index, Chem.MolToSmiles(mol_struct)]

    # if directory does not exist, create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(
        os.path.join(
            output_dir, "{}_{}.txt".format(id_smaller, id_bigger)
        ),
        "w",
    ) as f:
        # write space separated index and smiles in txt file
        for i, row in df.iterrows():
            f.write(str(row["atom_index"]) + " " + row["struct_smiles"] + "\n")


def main(args):
    disable_rdkit_logging()
    matches = get_data(args.matches_path, args.batch_number - 1, args.number_of_batches)
    # load the cached structures
    cached_structures = pickle.load(open(args.cached_structures_path, "rb"))
    # create the synthetic permutations
    for index, row in matches.iterrows():
        create_permutations(
            row["id_smaller"], row["id_bigger"], cached_structures, args.output_dir
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create cfmid synthetic permutations")
    parser.add_argument("matches_path", type=str)
    parser.add_argument("cached_structures_path", type=str)
    parser.add_argument("batch_number", type=int)
    parser.add_argument("number_of_batches", type=int)
    parser.add_argument("output_dir")

    args = parser.parse_args()
    main(args)

