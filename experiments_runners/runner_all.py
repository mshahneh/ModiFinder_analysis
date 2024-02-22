import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import pandas as pd
import argparse
import math
from runner_cfmid import main as run_cfmid
from runner_ours import main as run_ours
from runner_random import main as run_random
import json
import pickle

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

def run(id_known, id_modified, method, args, output, data_path):
    """
    Run the experiment
    """
    # call the proper runner based on the experiment method
    if method == "ours":
        run_ours(id_known, id_modified, method, args, output, data_path)
    elif method == "cfmid":
        run_cfmid(id_known, id_modified, method, args, output, data_path)
    elif method in ["random_choice", "random_distribution", "all_equal", "random_skewed", "multiple_random_choice", "multiple_random_distribution", "multiple_random_skewed"]:
        run_random(id_known, id_modified, method, args, output, data_path)
    else:
        raise ValueError("Invalid experiment method: {}".format(method))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('matches_path', type=str, help='path to matches.csv')
    parser.add_argument('experiments_arg', type=str, help='path to json file containing experiment arguments')
    parser.add_argument("data_path", type=str)
    parser.add_argument("batch_number", type=int)
    parser.add_argument("number_of_batches", type=int)
    parser.add_argument("output_dir")
    args = parser.parse_args()

    # load the experiment arguments
    experiment_args = json.load(open(args.experiments_arg, "r"))

    # load pairs
    matches = get_data(args.matches_path, args.batch_number - 1, args.number_of_batches)

    args = vars(args)
    args.update(experiment_args)

    # print(args.keys(), len(matches))
    for index, row in matches.iterrows():
        id_known = row["id_smaller"]
        id_modified = row["id_bigger"]
        output = os.path.join(args["output_dir"], "{}_{}.json".format(id_known, id_modified))
        run(id_known, id_modified, experiment_args["method"], args, output, args["data_path"])

# python runner_all.py  '/home/user/Substructure_Assignment/my_implementation/ModiFinder_Analysis/data/matches/GNPS-MSMLS.csv'  '/home/user/Substructure_Assignment/my_implementation/ModiFinder_Analysis/experiments_settings/0c7a48d0e14f845e6209983bcb98a65dc133780b.json' '/home/user/Substructure_Assignment/my_implementation/ModiFinder_Analysis/data' 0 1 '/home/user/Substructure_Assignment/my_implementation/ModiFinder_Analysis/results'