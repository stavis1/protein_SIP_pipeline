include { sipros } from './subworkflows/sipros'
include { isopacketModeler } from './subworkflows/isopacketModeler'
include { sipros_to_ipm } from './subworkflows/isopacketModeler'

params.results_dir = "$launchDir/results"

workflow {
    //make results directory
    file(params.results_dir).mkdir()
    file(params.design).copyTo(params.results_dir)

    Channel.of(file(params.design)).splitCsv(header : true, sep : '\t', strip : true)
        | sipros
        | sipros_to_ipm
        | isopacketModeler
}