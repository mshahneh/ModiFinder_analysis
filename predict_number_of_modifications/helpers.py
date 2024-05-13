import collections
import json
import requests
import numpy as np

SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["precursor_mz", "precursor_charge", "mz", "intensity"]
)

def convert_to_SpectrumTuple(peaks, precursor_mz, precursor_charge):
    """
    Converts the peaks to SpectrumTuple.
    params:
    peaks: list of peaks in the format [(mz, intensity), ...]
    precursor_mz: precursor m/z
    precursor_charge: precursor charge
    """
    if not isinstance(precursor_mz, (int, float)):
        try:
            precursor_mz = float(precursor_mz)
        except:
            raise ValueError('Precursor m/z must be a number')
    
    if not isinstance(precursor_charge, int):
        try:
            precursor_charge = int(precursor_charge)
        except:
            raise ValueError('Precursor charge must be an integer')
    
    peaks = normalize_peaks(peaks)
    res = {}
    res['precursor_charge'] = precursor_charge
    res['precursor_mz'] = precursor_mz
    res['mz'] = []
    res['intensity'] = []
    for peak in peaks:
        res['mz'].append(peak[0])
        res['intensity'].append(peak[1])
    
    return SpectrumTuple(**res)

def SpectrumTuple_to_dict(spectrum_tuple):
    """
    Converts the SpectrumTuple to a dictionary.
    """
    res = {}
    res['precursor_charge'] = spectrum_tuple.precursor_charge
    res['precursor_mz'] = spectrum_tuple.precursor_mz
    res['peaks'] = []
    for i in range(len(spectrum_tuple.mz)):
        res['peaks'].append((spectrum_tuple.mz[i], spectrum_tuple.intensity[i]))
    return res

def normalize_peaks(peaks):
    """
    l2 normalizes the peaks.
    """
    # l2 normalize the peaks over the intensity
    l2_norm = np.linalg.norm([peak[1] for peak in peaks])
    normalized_peaks = [(peak[0], peak[1] / l2_norm) for peak in peaks]
    return normalized_peaks

def get_usi_data(usi):
    """
    Requests the data using the USI.
    """
    url = 'https://metabolomics-usi.gnps2.org/json/' + "?usi1=" + usi
    try:
        r = requests.get(url, timeout=10)
    except requests.exceptions.Timeout:
        raise ValueError('Request timed out')
    
    data = json.loads(r.text)
    # change key names from precaursor_mz to Precursor_MZ
    data['Precursor_MZ'] = data.pop('precursor_mz')
    data['Charge'] = data.pop('precursor_charge')

    return data