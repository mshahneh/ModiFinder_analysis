import pandas as pd
import sys
import os
import argparse
import json
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from combine_csvs import main as combine_csvs_main
from eval import main as eval_main

def main(experiments_meta_path, skip_existing, just_eval=False):
    project_root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    run_config = json.load(open(os.path.join(project_root, "run_config.json")))
    data_folder = run_config["data_folder"]
    if not os.path.isabs(data_folder):
        data_folder = os.path.abspath(os.path.join(project_root, data_folder))

    results_dir = run_config["result_folder"]
    if not os.path.isabs(results_dir):
        results_dir = os.path.abspath(os.path.join(project_root, results_dir))
    
    df = pd.read_csv(experiments_meta_path)
    for i, row in df.iterrows():
        exp_id = row["id"]
        library_name = row["library"]
        method = row["method"]
        print(
            f"running experiment {i+1} out of {len(df)} with library: {library_name} and method: {method} fragment_filter: {row['fragment_filter']}"
        )

        matches_path = os.path.join(data_folder, "matches", library_name + ".csv")
        experiments_arg = os.path.join(project_root, "experiments_settings", exp_id + ".json")
        data_path = data_folder
        matches = pd.read_csv(matches_path)
        number_of_batches = min(100, len(matches)//3)


        directory_path = os.path.join(
            results_dir, exp_id, "evals"
        )
        output_path = os.path.join(
            results_dir, exp_id, "combined.csv"
        )

        if just_eval:
            # for all files in the directory, run eval
            cmnd2 = "nextflow run {}/experiments_runners/evals.nf  --matches_path {} --data_path {} --number_of_batches {} --results_dir {}".format(
                project_root,
                matches_path, os.path.join(results_dir, exp_id, "predictions"), number_of_batches, os.path.join(results_dir, exp_id)
            )
            os.system(cmnd2)
        else:
            if skip_existing and os.path.exists(output_path):
                print(f"skipping experiment {i+1} with library: {library_name}")
            else:
                if not os.path.exists(os.path.join(results_dir, exp_id)):
                    os.makedirs(os.path.join(results_dir, exp_id))
                cmnd = "nextflow run {}/experiments_runners/runner.nf  --matches_path {} --experiments_arg {} --data_path {} --number_of_batches {} --results_dir {}".format(
                    project_root,
                    matches_path, experiments_arg, data_path, number_of_batches, os.path.join(results_dir, exp_id)
                )
                # print(cmnd)
                os.system(cmnd)

                cmnd2 = "nextflow run {}/experiments_runners/evals.nf  --matches_path {} --data_path {} --number_of_batches {} --results_dir {}".format(
                    project_root,
                    matches_path, os.path.join(results_dir, exp_id, "predictions"), number_of_batches, os.path.join(results_dir, exp_id)
                )
                os.system(cmnd2)

        # combine all csv files in the directory

        combine_csvs_main(directory_path, output_path)
        print(f"done with experiment {i+1} with library: {library_name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run all experiments in a experiments meta file")
    parser.add_argument("experiments_meta_path",type=str)
    # add an option to skip existing experiments
    parser.add_argument("--skip_existing", action="store_true", default=False)
    parser.add_argument("--just_eval", action="store_true", default=False)
    args = parser.parse_args()
    main(args.experiments_meta_path, args.skip_existing, args.just_eval)
