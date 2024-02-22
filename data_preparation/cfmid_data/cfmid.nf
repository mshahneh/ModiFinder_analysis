#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// param to the data dir
params.data_dir = "/home/user/LabData/Reza/data"
params.library_name = "GNPS-MSMLS"
params.matches_path = "${params.data_dir}/matches/${params.library_name}.csv"
params.cached_structures_path = "${params.data_dir}/cached_structures/${params.library_name}.pkl"
params.batch_count = 10
params.cfmid_path = "/cfmid/public/"

process Sysntesis {
    errorStrategy 'ignore'
    maxForks 90
    publishDir "${params.data_dir}/cfmid_exp", mode: 'copy'

    input:
    val batch_id

    output:
    file "synth_output/*" optional true

    script:
    """
    mkdir synth_output
    python $baseDir/mol_struct_synth.py ${params.matches_path}  ${params.cached_structures_path}  ${batch_id}  ${params.batch_count} synth_output/
    """
}

process PerformCfmid{
    maxForks 90
    
    container 'wishartlab/cfmid:latest'
    shell = ['/bin/sh']
    publishDir "${params.data_dir}/cfmid_exp", mode: 'copy'
    // errorStrategy 'ignore'
    // time '5m'

    input:
    file input_synth

    output:
    file "cfmid_preds/*"

    """
    start=\${PWD}
    mkdir cfmid_preds
    cd $params.cfmid_path &&  cfm-predict \$start/$input_synth 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 \$start/cfmid_preds/$input_synth
    """
}

workflow {
    // create a channel for each batch
    def batches = Channel.from(1..params.batch_count)
    def synth_out = Sysntesis(batches).flatten()
    // synth_out.view()
    PerformCfmid(synth_out)
}