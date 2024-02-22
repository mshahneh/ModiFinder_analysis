#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.matches_path = ""
params.experiments_arg = ""
params.number_of_batches = ""
params.data_path = ""
params.results_dir = ""


process PerformEval {
    // for each file in the predictions folder, run the evaluation script and add the result as a line in a csv file
    // conda "${baseDir}/../environment.yml"
    conda "/home/user/mambaforge/envs/modi-finder-analysis"
    maxForks 90
    publishDir "${params.results_dir}", mode: 'copy'

    input:
    val batch_number

    output:
    file "evals/*.csv"

    script:
    """
    mkdir evals
    python $baseDir/eval.py ${params.matches_path} ${params.data_path} ${batch_number} ${params.number_of_batches} evals/
    """
}

workflow {
    def batches = Channel.from(1..params.number_of_batches)
    PerformEval(batches)
}