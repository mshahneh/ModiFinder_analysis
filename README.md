# ModiFinder Analysis

This repository provides codes and examples for the analysis of the paper:
``` ModiFinder: Tandem Mass Spectral Alignment Enables Structural Modification Site Localization ```

_Mohammad Reza Zare Shahneh, Michael Strobel, Giovanni Andrea Vitale, Christian Geibel, Vanessa V Phelan, Daniel Petras, Allegra T Aron, Yasin El Abiead, Neha Garg, Mingxun Wang_

## Install and setup
1. After cloning the repository, you need to add the ModiFinder module:

    ```git submodule update --init --recursive```

1. Install the conda enviroment, We recommend using mamba instead of conda for fast install (e.g., `mamba env create -f environment.yml`):

    ` conda env create -f environment.yml`

1. Install [`nextflow`](https://www.nextflow.io/docs/latest/getstarted.html)

1. Activate the environment:

    `conda activate modi-finder-analysis`


## Data
First, you need to set the directory of the data in `run_config` file. Then you can download the data or create it from scratch:
* You can download the files used in this project from: [`Zenodo`](https://zenodo.org/records/10674462?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjQwZWM4NjA0LTVmNDctNGFkNi1iNDgxLWZhNWE0NzMzYjBmMSIsImRhdGEiOnt9LCJyYW5kb20iOiJhNjE1ZjA1NGQ1MGY0MGQzNjk0MmU1YmZmNjg0NzAyMCJ9.NqwtedxTZyGK2Df2GgqU3Z2IMetDuSkFi5p7wprp0kzjHxse_w-KY3wlChw) and put them in the data directory defined earlier. The final format should be similar to this:
    ```
    your_data_directory/
    ├── matches/
    ├── helpers/
    ├── SIRIUS/
    └── cfmid_exp/
    ```

    **Please note that if you choose to download data in this manner, due to the necessity of requesting information for each individual compound in real-time, it is essential to restrict the number of concurrent processes to avoid exceeding the server's request limits.**

* You can download and create the data used in this project from scratch, by running the `data_prepare_main.py`:
    ```
    conda activate modi-finder-analysis
    python ./data_preparation/data_prepare_main.py
    ```

    **Please note that the data for SIRIUS has to be dowloaded from the provided link in the previous section or use `gnps2` to run the workflow.**

## Experiments
To run our experiments, you can run the following command:
    ```
    conda activate modi-finder-analysis
    python ./experiments_runners/experiments_runner.py './experiments_settings/all_experiments_settings.csv'
    ```

## Results
### Performance result
you can check `paper_figures/performance_results.ipynb` notebook for performance result illustrations.
### Helpers contribution
you can check  `paper_figures/how_much_helpers_help.ipynb` notebook for helpers contribution.
### Evaluation score illustration
you can check `paper_figures/evaluation_score_illustration.ipynb` notebook for evaluation score illustration.
### dataset stats
you cak check `paper_figures/datasets.ipynb` notebook for dataset stats.
