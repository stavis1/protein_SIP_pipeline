params.results_dir = "$launchDir/results"

process sipros_config_generator {
    container 'stavisvols/psp_sipros:latest'
    containerOptions "--bind $launchDir:/data/"
    publishDir params.results_dir, mode: 'copy'

    input:
    val row

    output:
    path '*.cfg'

    script:
    """
    conda run -n sipros_env python configGenerator -i $row.sipros_config -o ./ -e $row.label_elm
    """
}

process sipros_convert_raw_file {
    container 'stavisvols/psp_sipros:latest'
    containerOptions "--bind $launchDir:/data/"

    input:
    val row

    output:
    path '*.FT{1,2}'

    script:
    //figure out how you want to determine the number of allocated cores then pass that to Raxport with the -j flag
    """
    conda run -n sipros_env mono /opt/conda/envs/sipros_env/bin/Raxport.exe -i ./ -o ./
    """
}

process sipros_search {
    container 'stavisvols/psp_sipros:latest'
    containerOptions "--bind $launchDir:/data/"

    input:
    tuple path(congif_file), path(ft_files), val(row)

    output:
    path '*.sip'

    script:
    """
    conda run -n sipros_env 
    """
}

process sipros_to_ipm{

}

process ipm_parse_data{

}

workflow sipros {
    take:
    row

    main:
    //set up per-file data as value channels
    row_channel = channel.value(row)
    ft_files = sipros_convert_raw_file(row_channel)
        | collect

    //run searches at each % RIA step
    config_files = sipros_config_generator(row_channel)
        | flatten
    search_results = sipros_search(config_files, ft_files, row_channel)
        | collect
    
    //do post-processing
    processed_results = 

    emit:
    getLabelPCTinEachFT.out
}

workflow {
    //make results directory
    file(params.results_dir).mkdir()
    file(params.design).copyTo(params.results_dir)

    //parse the design file
    design = Channel.of(file(params.design)).splitCsv(header : true, sep : '\t', strip : true)

    //do per-file processing
    ipm_step_1 = sipros(design)
        | sipros_to_ipm
        | ipm_parse_data
        | collect
    
    //run IPM classifier step

    //run IPM fitting jobs

    //
}