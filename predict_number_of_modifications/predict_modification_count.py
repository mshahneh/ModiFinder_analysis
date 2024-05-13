import argparse
import json
import requests
import math
import numpy as np
import methods as methods
import pickle
import copy
from helpers import *

def handle_spectrum(spectrum):
    """
    Handles the spectrum input and returns the SpectrumTuple.
    """
    if isinstance(spectrum, SpectrumTuple):
        return copy.deepcopy(spectrum)
    elif isinstance(spectrum, dict):
        data = copy.deepcopy(spectrum)
    elif isinstance(spectrum, str):
        if spectrum.endswith('.mgf'):
            raise ValueError('MGF file format not supported')
        elif spectrum.endswith('.mzML'):
            raise ValueError('mzML file format not supported')
        elif spectrum.endswith('.mzXML'):
            raise ValueError('mzXML file format not supported')
        elif spectrum.endswith('.json'):
            with open(spectrum, 'r') as f:
                data = json.load(f)
        elif spectrum.endswith('.pkl'):
            with open(spectrum, 'rb') as f:
                data = pickle.load(f)
        else:
            data = get_usi_data(spectrum)
    else:
        raise ValueError('Unknown spectrum input format')
    
    if isinstance(data, SpectrumTuple):
        return data
    
    if 'peaks' not in data or 'Precursor_MZ' not in data or 'Charge' not in data:
        raise ValueError('Unknown spectrum input format')
    return convert_to_SpectrumTuple(data['peaks'], data['Precursor_MZ'], data['Charge'])

def main(spectrum1 , spectrum2, method, run_args):
    spectrum1_data = handle_spectrum(spectrum1)
    spectrum2_data = handle_spectrum(spectrum2)
    if "ppm" not in run_args:
        run_args["ppm"] = 40
    if "mz_tolerance" not in run_args:
        run_args["mz_tolerance"] = 0.1

    thresh = min(run_args['ppm'] * spectrum1_data.precursor_mz / 1e6, run_args['mz_tolerance'])
    if method == 'exhaustive':
        return methods.approximate_exhaustive(spectrum1_data, spectrum2_data, thresh)
    else:
        raise ValueError('Unknown method')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Predict the number of modifications in a spectrum')
    parser.add_argument('spectrum1', help='The first spectrum')
    parser.add_argument('spectrum2', help='The second spectrum')
    parser.add_argument('method', help='The method to use for prediction', default='exhaustive')
    parser.add_argument('--ppm', type=float, help='The ppm to use for prediction', default=40)
    parser.add_argument('--mz_tolerance', type=float, help='The m/z tolerance to use for prediction', default=0.01)

    args = parser.parse_args()

    run_args = {
        'ppm': args.ppm,
        'mz_tolerance': args.mz_tolerance
    }

    main(args.spectrum1, args.spectrum2, args.method, run_args)

# example:
# python predict_modification_count.py 'mzspec:GNPS:TUEBINGEN-NATURAL-PRODUCT-COLLECTION:accession:CCMSLIB00011906190' 'mzspec:GNPS:TUEBINGEN-NATURAL-PRODUCT-COLLECTION:accession:CCMSLIB00011906105' exhaustive