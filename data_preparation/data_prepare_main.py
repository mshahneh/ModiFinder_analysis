# scriptDir=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
# # download and prepare the data 
# python $scriptDir/prepare_datasets.py
# #find matches and helpers
# python $scriptDir/prepare_matches.py
# #calculate cfmid pairs
# python $scriptDir/calculate_cfmid_predictions.py

from urllib.request import urlopen
import pandas as pd
import os
import json

import prepare_datasets as prepare_datasets
import prepare_matches as prepare_matches
import calculate_cfmid_predictions as calculate_cfmid_predictions

project_root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
run_config = json.load(open(os.path.join(project_root, "run_config.json")))
data_folder = run_config["data_folder"]
# if data_folder is not absolute path, make it absolute
if not os.path.isabs(data_folder):
    data_folder = os.path.abspath(os.path.join(project_root, data_folder))

# # download the libraries
# prepare_datasets.main(project_root, run_config)

# # find matches and helpers
prepare_matches.main(project_root, run_config)

# combine all the matches
all_matches = pd.DataFrame()
libraries = pd.read_csv(os.path.join(project_root, "libraries.csv"))
for i, row in libraries.iterrows():
    library_compounds_data = pd.read_csv(os.path.join(data_folder, "libraries", row['Library'] + "_compounds_data.csv"))
    items_to_remove = library_compounds_data[library_compounds_data['num_atoms'] > 100]['spectrum_id']


    library = row["Library"]
    matches = pd.read_csv(os.path.join(data_folder, "matches", library + ".csv"))
    new_matches = matches[~matches['id_bigger'].isin(items_to_remove)]
    new_matches = new_matches[~new_matches['id_smaller'].isin(items_to_remove)]

    new_matches.to_csv(os.path.join(data_folder, "matches", library + ".csv"), index=False)
    # add column with library name
    new_matches["library"] = library
    all_matches = pd.concat([all_matches, new_matches])
all_matches.to_csv(os.path.join(data_folder, "matches", "all_matches.csv"), index=False)


# calculate cfmid pairs
# calculate_cfmid_predictions.main(project_root, run_config)