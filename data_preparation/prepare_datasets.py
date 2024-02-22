"""
This file downloads the libraries from the given links, cleans the data and caches it.
"""

from urllib.request import urlopen
import pandas as pd
import os
import json
from rdkit import Chem
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
from tqdm import tqdm
import pickle

def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')

def download(Library, Link, data_folder, skip_download):
    """
    download a library from a given link
    """
    
    print(f"Downloading {Library} ...", end=" ")
    
    if skip_download:
        try:
            with open(os.path.join(data_folder, "libraries", Library + ".json"), "r") as f:
                library_data = f.read()
            print("Done")
            return library_data
        except:
            print(f"Error handling {Library}.json . Skipping...")
            return None
    
    response = urlopen(Link)
    data_json = json.loads(response.read())

    if not os.path.exists(os.path.join(data_folder, "libraries")):
        os.makedirs(os.path.join(data_folder, "libraries"))

    # store the content
    json.dump(data_json, open(os.path.join(data_folder, "libraries", Library + ".json"), "w"))
    print("Done")
    return data_json

def extract_data(compound_data, columns):
    row_data = {}
    for col in columns:
        try:
            row_data[col] = compound_data[col]
        except:
            row_data[col] = None


def clean_and_cache(library_data, data_folder, library_name, cleaned_spectras):
    print(f"cleaning and caching {library_name} ...")
    count_not_clean = 0
    count_not_smiles = 0
    
    compounds_data = {}
    structures_data = {}

    columns = ["spectrum_id", "num_aromatic_rings", "num_atoms", "num_bonds", "num_rings", "not_connected", "formula_smiles"]
    structures_meta_data = pd.DataFrame(columns=columns)

    for i in tqdm(range(len(library_data)), ascii=True):
        compoundID = library_data[i]['spectrum_id']

        if (cleaned_spectras is not None) and (compoundID not in cleaned_spectras):
            count_not_clean += 1
            continue
        try:
            structure = Chem.MolFromSmiles(library_data[i]['Smiles'])
            # print(library_data[i]['Smiles'], structure)
            if structure is None:
                count_not_smiles += 1
                continue

            compound_data = {}
            compound_data['spectrum_id'] = compoundID
            compound_data['num_aromatic_rings'] = rdMolDescriptors.CalcNumAromaticRings(structure)
            compound_data['num_atoms'] = structure.GetNumAtoms()
            compound_data['num_bonds'] = structure.GetNumBonds()
            compound_data['num_rings'] = rdMolDescriptors.CalcNumRings(structure)
            compound_data['not_connected'] = 1 if "." in library_data[i]['Smiles'] else 0
            compound_data['Formula_smiles'] = library_data[i]['Formula_smiles']
            structures_meta_data = pd.concat([structures_meta_data, pd.DataFrame(compound_data, index=[0])], ignore_index=True)
        except:
            count_not_smiles += 1
            continue

        if not os.path.exists(os.path.join(data_folder, "cached_compounds")):
            os.makedirs(os.path.join(data_folder, "cached_compounds"))
        with open(os.path.join(data_folder, "cached_compounds", compoundID + ".pkl"), "wb") as f:
            pickle.dump(library_data[i], f)
        
        if not os.path.exists(os.path.join(data_folder, "cached_structures")):
            os.makedirs(os.path.join(data_folder, "cached_structures"))
        with open(os.path.join(data_folder, "cached_structures", compoundID + ".pkl"), "wb") as f:
            pickle.dump(structure, f)
        
        compounds_data[compoundID] = library_data[i]
        structures_data[compoundID] = structure

    print(f"Done. {count_not_clean} not clean, {count_not_smiles} not smiles")

    with open(os.path.join(data_folder, "cached_compounds", library_name + ".pkl"), "wb") as f:
        pickle.dump(compounds_data, f)
    
    with open(os.path.join(data_folder, "cached_structures", library_name + ".pkl"), "wb") as f:
        pickle.dump(structures_data, f)
    
    structures_meta_data.to_csv(os.path.join(data_folder, library_name + "_structures_meta.csv"), index=False)
    
    meta_data = {}
    meta_data["count_not_smiles"] = count_not_smiles
    meta_data["count_not_clean"] = count_not_clean
    meta_data["count_total"] = len(library_data)

    print("Done")
    return meta_data

def main(project_root, run_config):
    disable_rdkit_logging()
    meta_data = pd.DataFrame()
    
    data_folder = run_config["data_folder"]
    # if data_folder is not absolute path, make it absolute
    if not os.path.isabs(data_folder):
        data_folder = os.path.abspath(os.path.join(project_root, data_folder))

    libraries = pd.read_csv(os.path.join(project_root, "libraries.csv"))
    
    df = pd.read_csv(os.path.join(project_root, "ALL_GNPS_cleaned.csv"))
    df = df[df['ppmBetweenExpAndThMass'] < run_config['data_prepare']['clean_threshold_ppm']]
    cleaned_spectras = set(df["spectrum_id"].values)
    
    for i, row in libraries.iterrows():
        library_meta = {}
        library_meta["Library"] = row['Library']

        library_data = download(row['Library'], row['Link'], data_folder, run_config['data_prepare']['skip_download'])
        if not run_config['data_prepare']['skip_cache']:
            temp_meta = clean_and_cache(library_data, data_folder, row['Library'], cleaned_spectras)
            library_meta.update(temp_meta)
            meta_data = pd.concat([meta_data, pd.DataFrame(library_meta, index=[0])], ignore_index=True)
        
        # read structures_meta_data
        structures_meta_data = pd.read_csv(os.path.join(data_folder, row['Library'] + "_structures_meta.csv"))
        # combine with cleaned_data
        library_compounds_data = df.merge(structures_meta_data, on="spectrum_id", how="inner")
        # select only the columns we need
        columns = ['spectrum_id', 'collision_energy', 'Adduct', 'Compound_Source', 'Formula_smiles',
        'Compund_Name', 'Precursor_MZ', 'Precursor_MZ', 'Charge', 'msMassAnalyzer', "num_aromatic_rings", "num_atoms", "num_bonds", "num_rings", "not_connected"]
        library_compounds_data = library_compounds_data[columns]

        library_meta["mass_bigger_than_threshold"] = len(library_compounds_data[library_compounds_data['Precursor_MZ'] > run_config['data_prepare']['mass_threshold']])
        library_compounds_data = library_compounds_data[library_compounds_data['Precursor_MZ'] <= run_config['data_prepare']['mass_threshold']]

        if not run_config['data_prepare']['count_not_connected']:
            library_meta['count_not_connected'] = library_compounds_data['not_connected'].sum()
            library_compounds_data = library_compounds_data[library_compounds_data['not_connected'] == 0]

        # store the data
        library_compounds_data.to_csv(os.path.join(data_folder, "libraries", row['Library'] + "_compounds_data.csv"), index=False)

    meta_data.to_csv(os.path.join(data_folder, "meta_data.csv"), index=False)

if __name__ == "__main__":
    project_root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    run_config = json.load(open(os.path.join(project_root, "run_config.json")))
    main(project_root, run_config)
    




        
        


    


