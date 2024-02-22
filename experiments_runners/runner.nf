#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.matches_path = ""
params.experiments_arg = ""
params.number_of_batches = ""
params.data_path = ""
params.results_dir = ""


process CalcProcess {
    // errorStrategy 'ignore'
    // conda "${baseDir}/../environment.yml"
    conda "/home/user/mambaforge/envs/modi-finder-analysis"
    maxForks 90
    publishDir "${params.results_dir}", mode: 'copy'

    input:
    val batch_number
    output:
    file "predictions/*"

    script:
    """
    mkdir predictions
    python $baseDir/runner_all.py ${params.matches_path} ${params.experiments_arg} ${params.data_path} ${batch_number} ${params.number_of_batches} predictions/
    """
}

workflow {
    def batches = Channel.from(1..params.number_of_batches)
    def predictions = CalcProcess(batches)
}