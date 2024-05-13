from helpers import *
from typing import List, Tuple
import numpy as np

def make_graph(spec: SpectrumTuple,
                spec_other: SpectrumTuple,
                fragment_mz_tolerance: float,
                fragment_ppm_tolerance: float,
                list_diffs: List[float]
               ) -> np.ndarray:
    """
    Make a adjacency matrix for the graph.
    params:
    spec: SpectrumTuple object for the first spectrum
    spec_other: SpectrumTuple object for the second spectrum
    fragment_mz_tolerance: float value for the m/z tolerance in Da
    fragment_ppm_tolerance: float value for the ppm tolerance in ppm
    list_diffs: list of possible differences in m/z
    """
    # initialize the graph
    graph = np.zeros((len(spec.mz), len(spec_other.mz)))

    for peak_index, (peak_mz, peak_intensity) in enumerate(zip(spec.mz, spec.intensity)):
        for peak_other_index, (peak_other_mz, peak_other_intensity) in enumerate(zip(spec_other.mz, spec_other.intensity)):
            peak_diff = peak_mz - peak_other_mz
            for diff in list_diffs:
                if (abs(peak_diff - diff) <= fragment_mz_tolerance) and ((abs(peak_diff - diff) / peak_mz) * 1e6 <= fragment_ppm_tolerance):
                    graph[peak_index][peak_other_index] = max(graph[peak_index][peak_other_index], peak_intensity * peak_other_intensity)
                    break
    return graph


def greedy_matching(graph: np.ndarray) -> List[Tuple[int, int]]:
    """
    Calculate the matching that maximizes the cosine similarity using the greedy method.
    params:
    graph: adjacency matrix for the graph
    """
    score, peak_matches = 0.0, []
    values = []
    for i in range(graph.shape[0]):
        for j in range(graph.shape[1]):
            if graph[i][j] > 0:
                values.append((graph[i][j], (i, j)))
    values = sorted(values, key=lambda x: x[0], reverse=True)
    used_rows, used_columns = np.zeros(graph.shape[0]), np.zeros(graph.shape[1])
    for value, (i, j) in values:
        if used_rows[i] == 0 and used_columns[j] == 0:
            score += value
            peak_matches.append((i, j))
            used_rows[i] = 1
            used_columns[j] = 1
    return score, peak_matches


def calculate_matching(graph: np.ndarray, method: str) -> Tuple[float, List[Tuple[int, int]]]:
    """
    Calculate the matching that maximizes the cosine similarity.
    params:
    graph: adjacency matrix for the graph
    method: method to use for matching
    """
    if method == 'greedy':
        return greedy_matching(graph)
    elif method == 'hungarian':
        pass
    else:
        raise ValueError('Invalid method')

def cosine_slow(
    spec: SpectrumTuple,
    spec_other: SpectrumTuple,
    fragment_mz_tolerance: float,
    fragment_ppm_tolerance: float,
    allow_shift: bool = False,
    list_diffs: List[float] = None,
) -> Tuple[float, List[Tuple[int, int]]]:
    """
    Calculate the cosine similarity between two spectra.
    params:
    spec: SpectrumTuple object for the first spectrum
    spec_other: SpectrumTuple object for the second spectrum
    fragment_mz_tolerance: float value for the m/z tolerance in Da
    fragment_ppm_tolerance: float value for the ppm tolerance in ppm
    allow_shift: when list_diffs is none, this parameter is used to allow matching with precursor difference shifts
    list_diffs: list of possible differences in m/z
    """
    if spec.precursor_charge != 1 or spec_other.precursor_charge != 1:
        raise ValueError("Only precursor charge 1 is supported")
    
    if list_diffs is None:
        list_diffs = list([0])
        if allow_shift:
            list_diffs.append(spec_other.precursor_mz - spec.precursor_mz)

    graph = make_graph(spec, spec_other, fragment_mz_tolerance, fragment_ppm_tolerance, list_diffs)

    return calculate_matching(graph, method='greedy')