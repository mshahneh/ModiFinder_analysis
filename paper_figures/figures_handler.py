# do the imports
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
import seaborn as sns
from scipy.stats import gaussian_kde
from textwrap import wrap

# add parent directory to path
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from SmallMol_Mod_Site_Localization import Compound_n as Compound
from SmallMol_Mod_Site_Localization import ModificationSiteLocator as modSite
from SmallMol_Mod_Site_Localization import utils_n as utils
from SmallMol_Mod_Site_Localization import handle_network as handle_network
from SmallMol_Mod_Site_Localization import visualizer
from SmallMol_Mod_Site_Localization import calculate_scores_n as calculate_scores
from SmallMol_Mod_Site_Localization import alignment_n as alignment


def read_library(library_name, path):
    df = pd.read_csv(path)
    df = df[df['library'] == library_name]
    return df

def read_data_of_runs(path, df):
    data = {}
    for i, row in df.iterrows():
        id = row['id']
        temp_df = pd.read_csv(os.path.join(path, id,  "combined.csv"))
        data[id] = temp_df
    
    return data

def get_method_to_id_map(df):
    methods = {}
    for i, row in df.iterrows():
        # if fragment_filter is not defined or na, skip
        if row["fragment_filter"] == "None" or row["fragment_filter"] == None or row["fragment_filter"] == "NaN" or row["fragment_filter"] == "nan" or pd.isna(row["fragment_filter"]):
            methods[row["method"]] = row['id']
        else:
            methods[row["fragment_filter"]] = row['id']
            # print("processing", row["id"], i, row["fragment_filter"], pd.isna(row["fragment_filter"]))

    return methods

def get_all(library_name, path_to_experiment, path_to_data):
    df = read_library(library_name, path_to_experiment)
    data = read_data_of_runs(path_to_data, df)
    methods = get_method_to_id_map(df)
    return df, data, methods

def merge_dataframes(res, df2, name, cols_compare, cols_val, cols_meta):
    cols_to_keep = []
    for col in cols_meta:
        if col in df2.columns and col not in res.columns:
            cols_to_keep.append(col)
    
    # rename cols in cols_val with method name as prefix
    if isinstance(cols_val, str): # if cols_compare is instance of str, then it is just one column and just use name
        # if there is a column with the same name in df2, remove it
        if name in df2.columns:
            df2 = df2.drop(columns=[name])
        cols_val_renamed = name
        df2 = df2.rename(columns={cols_val: cols_val_renamed})
        cols_val_renamed = [cols_val_renamed]
    else:
        cols_val_renamed = []
        for col in cols_val:
            cols_val_renamed.append("{}_{}".format(name, col))
        df2 = df2.rename(columns=dict(zip(cols_val, cols_val_renamed)))
    
    # merge the dataframes
    merged = res.merge(df2[cols_compare + cols_val_renamed + cols_to_keep], on=cols_compare, how="inner")
    return merged

def create_merged_libraries(data, methods, cols_compare, cols_val, cols_meta):
    # get a random method to start with
    if "oracke" in methods:
        random_method = "oracle"
    else:
        random_method = list(methods.keys())[0]
    res = data[methods[random_method]][cols_compare]
    for method in methods:
        res = merge_dataframes(res, data[methods[method]], method, cols_compare, cols_val, cols_meta)
    
    return res

def combine_merged_libraries(library_names, path_to_experiment, path_to_data):
    methods_data = {}
    for library in library_names:
        df, data, methods = get_all(library, path_to_experiment, path_to_data)
        for method in methods:
            if method not in methods_data:
                methods_data[method] = data[methods[method]]
            else:
                methods_data[method] = pd.concat([methods_data[method], data[methods[method]]])
    return methods_data

def draw_violin(cols_to_draw, cols_name, result_df, figure_name, font_size=18, ax = None):
    if ax is None:
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111)
    # show gray grid
    # ax.grid(color='#cccccc', linestyle='-', linewidth=0.1)

    for index, method in enumerate(cols_to_draw):
        # drop na values
        temp = result_df.dropna(subset=[method])
        values = temp[method]
        # draw violin plot for each method and show mean value
        ax.violinplot(values, positions=[index], showmeans=True, showextrema=False)
        ax.text(index, np.mean(values), round(np.mean(values), 3), ha='center', va='bottom', fontsize=font_size)

    # set y axis limit
    if max(result_df[cols_to_draw].max()) <= 1:
        ax.set_ylim([0, 1])

    if cols_name is not None:
        ax.set_xticks(range(len(cols_to_draw)))
        ax.set_xticklabels(cols_name, rotation=45, ha='right', fontsize=font_size)
    else:
        ax.set_xticks([])
        ax.set_xticklabels([])
    ax.tick_params(axis='y', labelsize=font_size)
    ax.set_title(figure_name, fontsize=font_size)

def draw_scatter(x_vals, y_vals, x_label, y_label, name):

    x = x_vals
    y = y_vals

    # # find indexes where x is 0
    # idx = np.where(y<0.1)[0]
    # # remove those indexes
    # x = np.delete(x, idx)
    # y = np.delete(y, idx)

    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    fig, ax = plt.subplots()
    fig.tight_layout()
    title = ax.set_title("\n".join(wrap("Evaluation Scores for {} where there is at least one shifted peak".format(name), 40)))
    title.set_y(1.05)
    fig.subplots_adjust(top=0.8)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    density = ax.scatter(x, y, c=z, s=10)

    # set y axis limit
    ax.set_ylim([0, 1])
    ax.set_xlim([0, 1])
    # drawy y=x line
    ax.plot([0, 1], [0, 1], color="red", linestyle="--")
    # # draw linear regression line
    # ax.plot(x, np.poly1d(np.polyfit(x, y, 1))(x), color="blue", linestyle="--")

    # add colorbar
    fig.colorbar(density)

def draw_heatmap(x_vals, y_vals, x_label, y_label, name, ax= None):
    # select vals of msbuddy
    x = x_vals
    y = y_vals

    bins = np.linspace(0.0, 1.0, 6)
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    if (ax == None):
        fig, ax = plt.subplots()
        fig.tight_layout()
    title = ax.set_title("\n".join(wrap("Evaluation Scores for {} where there is at least one shifted peak".format(name), 70)))
    title.set_y(1.05)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    # show color legend
    im = ax.imshow(heatmap.T, extent=extent, origin='lower', aspect='auto', cmap=plt.cm.jet)
    fig.colorbar(im, ax=ax)
    # ax.imshow(heatmap.T, extent=extent, origin='lower', aspect='auto', cmap='viridis')

def render_svg_to_png(svg_content):
    """Render SVG content to PNG using cairosvg."""
    s_png = cairosvg.svg2png(bytestring=svg_content)
    s_img = mpimg.imread(io.BytesIO(s_png))
    return s_img

def draw_the_peaks(ax, main_compound, mod_compound, modSiteLocator):
    # create the values for the bars
    unmatched = {'x':[], 'y':[]}
    matched_shifted = {'x':[], 'y':[]}
    matched_unshifted = {'x':[], 'y':[]}
    for peak in main_compound.peaks:
        unmatched['x'].append(peak[0])
        unmatched['y'].append(peak[1])
    for peak in mod_compound.peaks:
        unmatched['x'].append(peak[0])
        unmatched['y'].append(-peak[1])

    for shift in modSiteLocator.shifted:
        matched_shifted['x'].append(main_compound.peaks[shift[0]][0])
        matched_shifted['y'].append(main_compound.peaks[shift[0]][1])
        matched_shifted['x'].append(mod_compound.peaks[shift[1]][0])
        matched_shifted['y'].append(-mod_compound.peaks[shift[1]][1])
    for match in modSiteLocator.unshifted:
        matched_unshifted['x'].append(main_compound.peaks[match[0]][0])
        matched_unshifted['y'].append(main_compound.peaks[match[0]][1])
        matched_unshifted['x'].append(mod_compound.peaks[match[1]][0])
        matched_unshifted['y'].append(-mod_compound.peaks[match[1]][1])

    # draw the plot
    ax.set_title("peaks")
    ax.bar(unmatched['x'], unmatched['y'], color="grey")
    ax.bar(matched_shifted['x'], matched_shifted['y'], color="red")
    ax.bar(matched_unshifted['x'], matched_unshifted['y'], color="blue")
    # draw horizontal line at 0
    ax.axhline(0, color='black')
    # draw a vertical black dashed line at precursor mass that starts from 0 and ends at the highest peak
    ax.plot([main_compound.Precursor_MZ, main_compound.Precursor_MZ], [0, max(unmatched['y'])], color='black', linestyle='dashed')
    ax.plot([mod_compound.Precursor_MZ, mod_compound.Precursor_MZ], [0, min(unmatched['y'])], color='black', linestyle='dashed')

# draw the pairs
def draw_molecules(main_compound, mod_compound, modifLoc, fig = None):
    if fig is None:
        fig, ax = plt.subplots(1, 2, figsize=(10, 10))
    else:
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
        ax = [ax1, ax2]
    ax[0].set_title("modified compound")
    svg = visualizer.molToSVG(mod_compound.structure, main_compound.structure, True)
    ax[0].imshow(render_svg_to_png(svg))
    ax[1].set_title("main compound")
    svg = visualizer.highlightMolIndices(main_compound.structure, modifLoc)
    ax[1].imshow(render_svg_to_png(svg))

#draw fragments of a peak
def draw_frags_of_peak(modSiteLocator, peak_index, fig = None):
    print(modSiteLocator.get_structures_by_peak_index(peak_index))
    structures, locations, frags = modSiteLocator.get_structures_by_peak_index(peak_index)
    subplotsRows = math.ceil(len(locations)/3)
    if fig is None:
        fig, ax = plt.subplots(subplotsRows, 3, figsize=(10, 5))
    else:
        ax = []
        for i in range(subplotsRows):
            row = []
            for j in range(3):
                row.append(fig.add_subplot(subplotsRows, 3, i*3+j+1))
            ax.append(row)
        ax = np.array(ax)
        if subplotsRows == 1:
            ax = ax[0]
    for i, location in enumerate(locations):
        svgText = visualizer.highlightMolIndices(modSiteLocator.main_compound.structure, location)
        img = render_svg_to_png(svgText)
        if subplotsRows == 1:
            ax[i].imshow(img)
            # set title to fragment
            ax[i].set_title("fragment " + str(frags[i]))

        else:
            ax[i//3, i%3].imshow(img)
            # set title to fragment
            ax[i//3, i%3].set_title("fragment " + str(frags[i]))

def create_locator(id_smaller, id_bigger, add_adduct = True, args = {}):
    main_compound = compound.Compound(id_smaller, args=args)
    # try:
    #     with open(os.path.join(data_folder, "SIRIUS",  "{}_fragmentationtree.json".format(id_smaller)), "rb") as f:
    #             sirius = json.load(f)
    #     main_compound.apply_sirius(sirius, add_adduct)
    # except:
    #      print("no sirius data for", id_smaller)
    #      pass
    mod_compound = compound.Compound(id_bigger, args=args)
    modSiteLocator = modSite.ModificationSiteLocator(main_compound, mod_compound, args)
    true_modification_site = utils.calculateModificationSites(mod_compound.structure, main_compound.structure, False)[0]
    return modSiteLocator, true_modification_site, main_compound, mod_compound

def draw_match_info(modSiteLocator, true_modification_site, add_adduct = True):

    #draw 3 figures under each other
    fig1 = plt.figure(figsize=(10, 5))
    draw_molecules(modSiteLocator.main_compound, modSiteLocator.modified_compound, [true_modification_site], fig1)
    scores = modSiteLocator.generate_probabilities()
    evals = modSiteLocator.calculate_score(true_modification_site, "average_dist_normalized", scores)
    fig2, ax2 = plt.subplots(1, 2, figsize=(10, 3))
    ax2[0].set_title("scores _ with unshifted" + str(round(evals,2)))
    ax2[0].imshow(render_svg_to_png(visualizer.highlightScores(modSiteLocator.main_compound.structure, scores)))
    scores2 = modSiteLocator.generate_probabilities(shifted_only=True)
    evals2 = modSiteLocator.calculate_score(true_modification_site, "average_dist_normalized", scores2)
    ax2[1].set_title("scores - shifted only" + str(round(evals2,2)))
    ax2[1].imshow(render_svg_to_png(visualizer.highlightScores(modSiteLocator.main_compound.structure, scores2)))

    # print(modSiteLocator.shifted)
    # print(main_compound.Precursor_MZ, mod_compound.Precursor_MZ)
    for peak in modSiteLocator.shifted:
        print(modSiteLocator.modified_compound.peaks[peak[1]], modSiteLocator.main_compound.peaks[peak[0]], modSiteLocator.main_compound.peak_fragments_map[peak[0]])
    # print(modSiteLocator.cosine)
    fig3, ax3 = plt.subplots(1, 1, figsize=(10, 3))
    draw_the_peaks(ax3, modSiteLocator.main_compound, modSiteLocator.modified_compound, modSiteLocator)

def print_sirius_data(id_smaller, peak_weight):
    with open(os.path.join(data_folder, "SIRIUS",  "{}_fragmentationtree.json".format(id_smaller)), "rb") as f:
            sirius = json.load(f)
    print("peak_weight", peak_weight)
    index = utils.find_mz_in_sirius(sirius["fragments"], float(peak_weight), 0.01, 40)
    return sirius["fragments"][index]["molecularFormula"]

def print_match_info(ModiFinder, true_modification_site):
    temp_data = {}
    smaller_id = ModiFinder.main_compound.accession
    bigger_id = ModiFinder.modified_compound.accession
    temp_data["id_base"] = smaller_id
    temp_data["id_analog"] = bigger_id
    temp_data["correctly_predicted"] = (ModiFinder.calculate_score(true_modification_site, "is_max") > 0)
    temp_data["evaluation_score"] = ModiFinder.calculate_score(true_modification_site, "average_dist_normalized")
    temp_data["base_name"] = ModiFinder.main_compound.name
    temp_data["analog_name"] = ModiFinder.modified_compound.name
    temp_data["base_smiles"] = ModiFinder.main_compound.Smiles
    temp_data["analog_smiles"] = ModiFinder.modified_compound.Smiles
    temp_data["num_shifted_peaks"] = len(ModiFinder.shifted)
    temp_data["num_matches"] = len(ModiFinder.matched_peaks)
    temp_data["num_peaks_base"] = len(ModiFinder.main_compound.peaks)
    temp_data["delta_mass_mz"] = abs(ModiFinder.main_compound.Precursor_MZ - ModiFinder.modified_compound.Precursor_MZ)
    temp_data["link"] = handle_network.create_link_from_accession(smaller_id, bigger_id)

    # print data in 2 columns
    for key in temp_data:
        # align the keys and values
        print(key.ljust(30), temp_data[key])

def apply_helpers(ModiFinder, data_path = "/home/user/LabData/Reza/data"):
    id_known = ModiFinder.main_compound.accession
    id_modified = ModiFinder.modified_compound.accession
    with open(os.path.join(data_path, "helpers", "{}_{}.json".format(id_known, id_modified)), "rb") as f:
        helpers_data = json.load(f)
        ModiFinder.apply_helpers(helpers_data, "intersection")

def print_frags_of_shifted_peaks(ModiFinder):
    for peak in ModiFinder.shifted:
        print(ModiFinder.modified_compound.peaks[peak[1]], ModiFinder.main_compound.peaks[peak[0]], ModiFinder.main_compound.peak_fragments_map[peak[0]])

def print_helpers(ModiFinder, data_path = "/home/user/LabData/Reza/data"):
    id_known = ModiFinder.main_compound.accession
    id_modified = ModiFinder.modified_compound.accession
    with open(os.path.join(data_path, "helpers", "{}_{}.json".format(id_known, id_modified)), "rb") as f:
        helpers_data = json.load(f)
    
    print(len(helpers_data), helpers_data)

    if len(helpers_data) > 0:
        fig, ax = plt.subplots(1, len(helpers_data), figsize=(len(helpers_data) * 2, 5))
        for index, helper in enumerate(helpers_data):
            c = compound.Compound(helper)
            svg = visualizer.molToSVG( ModiFinder.main_compound.structure, c.structure, True)
            if len(helpers_data) == 1:
                ax.imshow(render_svg_to_png(svg))
            else:
                ax[index].imshow(render_svg_to_png(svg))


colors = ['#f2bdbe', '#bbd5e8', '#bfe2bf', '#ffd8b6', '#d4c4fb', '#f9f9f9']
method_names = {
    "none": "MF-N",
    "helpers": "MF-R",
    "oracle": "MF-O",
    "sirius": "MF-S",
    "msbuddy": "MF-B",
    "cfmid": "CFM-ID",
    "multiple_random_choice": "RC",
    "multiple_random_distribution": "RD"
}

method_colors = {
    "none": colors[0],
    "helpers": colors[1],
    "oracle": colors[2],
    "cfmid": colors[3],
    "multiple_random_choice": colors[4],
    "multiple_random_distribution": colors[5]
}

def get_basic_data(project_root):
    # read the config file
    with open(os.path.join(project_root, 'run_config.json')) as f:
        config = json.load(f)

    # set up the directories
    data_folder = config['data_folder']
    if not os.path.isabs(data_folder):
        data_folder = os.path.join(project_root, data_folder)
        data_folder = os.path.abspath(data_folder)
    
    results_directory = config['result_folder']
    if not os.path.isabs(results_directory):
        results_directory = os.path.join(project_root, results_directory)
        results_directory = os.path.abspath(results_directory)

    matches_directory = os.path.join(data_folder, 'matches')
    
    # read the list of libraries
    libraries = pd.read_csv(os.path.join(project_root, "libraries.csv"))

    # create a map from library to short name
    library_names = {row['Library']: row['short_name'] for index, row in libraries.iterrows()}

    return data_folder, results_directory, matches_directory, libraries, library_names