from helpers import *
from alignment import cosine_slow

def calculate_all_differences(spec1: SpectrumTuple,
                                spec2: SpectrumTuple,
                                thresh: float):
    """
    Calculate all the differences between the peaks of two spectra. returns a list of differences and their number of occurences.
    params:
    spec1: SpectrumTuple
    spec2: SpectrumTuple
    thresh: float
    """
    diffrences = [[0, 0]]
    diff_precursors = spec1.precursor_mz - spec2.precursor_mz

    for main_c_peak in spec1.mz:
        for mod_c_peak in spec2.mz:
            
            diff_mz = main_c_peak - mod_c_peak

            # if the signs are different, skip (only find shifts in the same direction)
            if diff_mz * diff_precursors < 0 and abs(diff_mz) > thresh: 
                continue

            # if the difference is too small, but not small enough to be unshifted match, skip
            if abs(diff_mz) < 10 and abs(diff_mz) > thresh: 
                # print(diff_mz)
                continue

            # if the difference is too large, skip
            if abs(diff_mz) > abs(diff_precursors) + thresh:
                continue
            
            found = False
            for diffrance in diffrences:
                if abs(diff_mz - diffrance[0]) < thresh:
                    found = True
                    diffrance[1] += 1
                    break
            if not found:
                diffrences.append([diff_mz, 1])

    # sort the peaks based on the mz difference
    diffrences = sorted(diffrences, key=lambda x: abs(x[0]))
    return diffrences

def calculate_possibilities(diffrences: list, thresh: float, target: float, count: int):
    """
    Calculate the possibilities of a set of differences to sum up to a target value. the set must contain count number of differences. returns a list of possibilities.
    params:
    diffrences: a list of differences as [difference, count]
    thresh: float threshold to consider two values equal
    target: float target value to sum up to
    count: int number of differences to sum up to the target value
    """

    res = []
    if count == 0:
        return res
    elif count == 1:
        for difference in diffrences:
            if abs(difference[0] - target) < thresh:
                res.append(difference)
        return res
    
    elif count == 2:
        left = 0
        right = len(diffrences) - 1
        while left < right:
            if abs(diffrences[left][0]) < thresh: # if the smaller peak is 0 (too small) don't consider it
                left += 1
                continue

            if abs(diffrences[left][0] + diffrences[right][0] - target) < thresh:
                res.append([diffrences[left], diffrences[right]])
            if abs(diffrences[left][0] + diffrences[right][0]) < abs(target):
                left += 1
            else:
                right -= 1
        
        return res
    
    else:
        for i in range(len(diffrences)):
            copied = diffrences.copy()
            copied.pop(i)
            possibilities = calculate_possibilities(copied, thresh, target - diffrences[i][0], count - 1)
            if len(possibilities) > 0:
                for possibility in possibilities:
                    res.append([diffrences[i]] + possibility)
        return res

def approximate_exhaustive(tuple1: SpectrumTuple, tuple2: SpectrumTuple, thresh: float):
    """
    returns a boolean value if the two spectra have one modification site or not.
    params:
    tuple1: SpectrumTuple
    tuple2: SpectrumTuple
    thresh: float
    """
    
    diffrences = calculate_all_differences(tuple1, tuple2, thresh)
    possibilities = [None, None, None]
    possibilities[0] = calculate_possibilities(diffrences, thresh, 0, 1)
    possibilities[1] = calculate_possibilities(diffrences, thresh, tuple1.precursor_mz - tuple2.precursor_mz, 1)
    possibilities[2] = calculate_possibilities(diffrences, thresh, tuple1.precursor_mz - tuple2.precursor_mz, 2)
    maxVal = [0, 0, 1, 0, 0, 0]
    for j in range(1, len(possibilities)):
        for item in possibilities[j]:
            if isinstance(item[0], list):
                diffs = [_[0] for _ in item]
                occurences = sum([_[1] for _ in item])
            else:
                diffs = [item[0]]
                occurences = item[1]
            
            diffs += [0, tuple1.precursor_mz - tuple2.precursor_mz]

            score, matches = cosine_slow(tuple1, tuple2, 0.01, 40, list_diffs=diffs)
            coef = (len(matches)/maxVal[2] + 19) / 20
            if score <= 0.4:
                continue
            if score >= maxVal[1] * coef and occurences >=  maxVal[3]:
                maxVal = [j, score, len(matches), occurences, item, matches]
    
    if maxVal[0] == 1:
        return True
    else:
        return False