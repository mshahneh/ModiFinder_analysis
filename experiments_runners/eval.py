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

from SmallMol_Mod_Site_Localization import calculate_scores_n as calculate_scores


def get_data(matches_path, batch_number, number_of_batches):
    """
    Get the data for the current batch
    """
    matches = pd.read_csv(matches_path)
    batch_size = math.floor(len(matches) / number_of_batches)
    if batch_number == number_of_batches - 1:
        batch = matches.iloc[batch_number * batch_size :]
    else:
        batch = matches.iloc[batch_number * batch_size : (batch_number + 1) * batch_size]
    return batch

def main(prediction_dict, output_path, method=None):
    if method is None:
        method = [
            "is_max",
            "proximity",
            "average_dist_normalized",
            "sorted_rank"
        ]
    if isinstance(method, str):
        method = [method]

    scores = {}

    file_name = output_path.split("/")[-1].split(".")[0]

    if "_" in file_name:
        id_known = file_name.split("_")[0]
        id_unknown = file_name.split("_")[1]
        scores["id_known"] = id_known
        scores["id_unknown"] = id_unknown

    for m in method:
        scores[m] = calculate_scores.calculate(
            prediction_dict["graph"],
            prediction_dict["probabilities"],
            prediction_dict["true_site"],
            m,
        )
    
    # add metadata to scores
    # any key that is not ["graph", "probabilities", "true_site"] is considered metadata
    for key in prediction_dict.keys():
        if key not in ["graph", "probabilities", "true_site"]:
            scores[key] = prediction_dict[key]
            # if numpy array, convert to 1 line string
            if isinstance(scores[key], np.ndarray):
                scores[key] = str(scores[key])
            
            if isinstance(scores[key], list):
                scores[key] = str(scores[key])
            
            if isinstance(scores[key], str):
                scores[key] = scores[key].replace("\n", " ")
                scores[key] = scores[key].replace("\r", " ")

    if id_known == "CCMSLIB00010122119" and id_unknown == "CCMSLIB00010123721":
        print(scores)

    # print(scores)
    # write scores to csv file with just one row
    df = pd.DataFrame([scores])
    # dont write the header
    # df.to_csv(output_path, index=False, header=False)
    df.to_csv(output_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="evaluate the graph prediction based on scoring"
    )
    # parser.add_argument(
    #     "prediction_path", type=str, help="path to the prediction file"
    # )
    # parser.add_argument(
    #     "--output_path",
    #     type=str,
    #     default="scores.csv",
    #     help="path to the output file",
    # )

    parser.add_argument('matches_path', type=str, help='path to matches.csv')
    parser.add_argument("data_path", type=str)
    parser.add_argument("batch_number", type=int)
    parser.add_argument("number_of_batches", type=int)
    parser.add_argument("output_path")
    args = parser.parse_args()

    matches = get_data(args.matches_path, args.batch_number - 1, args.number_of_batches)

    for index, row in matches.iterrows():
        try:
            id_known = row["id_smaller"]
            id_modified = row["id_bigger"]
            

            with open(os.path.join(args.data_path, "{}_{}.json".format(id_known, id_modified)), "rb") as f:
                prediction = json.load(f)

            # json to dict
            prediction_dict = {}
            for key in prediction.keys():
                # if type is list, convert to numpy array
                if isinstance(prediction[key], list):
                    prediction_dict[key] = np.array(prediction[key])
                else:
                    prediction_dict[key] = prediction[key]

            # if id_known == "CCMSLIB00010122119":
            #     print(prediction_dict)

            main(prediction_dict, os.path.join(args.output_path, "{}_{}.csv".format(id_known, id_modified)))
        except:
            pass

# python eval.py --prediction_path "/home/user/Substructure_Assignment/my_implementation/SmallMol_Mod_Site_Localization_Analysis/results/d735eede6567b2832a2ca6ed0849394557b8dccc/predictions/CCMSLIB00005463552_CCMSLIB00005464419.json" --output_path "/home/user/Substructure_Assignment/my_implementation/SmallMol_Mod_Site_Localization_Analysis/results/scores.csv"
# python eval.py --prediction_path "/home/user/Substructure_Assignment/my_implementation/SmallMol_Mod_Site_Localization_Analysis/results/test.json" --output_path "/home/user/Substructure_Assignment/my_implementation/SmallMol_Mod_Site_Localization_Analysis/results/scores.csv"
