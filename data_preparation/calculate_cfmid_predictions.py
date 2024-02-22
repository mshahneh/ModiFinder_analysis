from urllib.request import urlopen
import pandas as pd
import os
import json
import sys

def main(project_root, run_config):
    data_folder = run_config["data_folder"]
    # if data_folder is not absolute path, make it absolute
    if not os.path.isabs(data_folder):
        data_folder = os.path.abspath(os.path.join(project_root, data_folder))
    
    cfmid_path = run_config["cfmid_path"]
    # if cfmid_path is not absolute path, make it absolute
    if not os.path.isabs(cfmid_path):
        cfmid_path = os.path.abspath(os.path.join(project_root, cfmid_path))

    libraries = pd.read_csv(os.path.join(project_root, "libraries.csv"))

    for i, row in libraries.iterrows():
        try:
            matches_path = os.path.join(data_folder, "matches", row["Library"] + ".csv")
            cached_structures_path = os.path.join(data_folder, "cached_structures", row["Library"] + ".pkl")
            batch_count = run_config["nf_batch_count"]

            matches = pd.read_csv(matches_path)
            if len(matches) < batch_count:
                batch_count = len(matches)
            
            print("Running cfmid for library: ", row["Library"], " with ", batch_count, " batches")
        except:
            print("Error: ", sys.exc_info()[0])
            continue

        cfmid_nf_path = os.path.join(project_root, "data_preparation", "cfmid_data", "cfmid.nf")
        cmnd = "nextflow run {} --data_dir {} "\
            "--library_name {} --matches_path {} --cached_structures_path {} "\
                "--batch_count {} --cfmid_path {}".format(cfmid_nf_path, 
                                          data_folder, row["Library"], 
                                          matches_path, cached_structures_path,
                                            batch_count, cfmid_path)
        os.system(cmnd)

if __name__ == "__main__":
    project_root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    run_config = json.load(open(os.path.join(project_root, "run_config.json")))
    main(project_root, run_config)