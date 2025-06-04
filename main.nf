params.results_dir = "$launchDir/results"
//required columns in design file : sample_ID, label_elm, raw_file, sipros_config

process sipros_config_generator {
    container 'stavisvols/psp_sipros:latest'
    label 'sipros_small'
    // publishDir params.results_dir, mode: 'copy'

    input:
    tuple val(row), path(config_file)

    output:
    path '*.cfg'

    script:
    """
    conda run -n sipros_env python configGenerator -i $config_file -o ./ -e $row.label_elm
    """
}

process sipros_convert_raw_file {
    container 'stavisvols/psp_sipros:latest'
    label 'sipros_med'

    input:
    path rawfile

    output:
    path '*.FT{1,2}'

    script:
    //figure out how you want to determine the number of allocated cores then pass that to Raxport with the -j flag
    """
    conda run -n sipros_env mono /opt/conda/envs/sipros_env/bin/Raxport.exe -i ./ -o ./ -j 3
    """
}

process sipros_search {
    container 'stavisvols/psp_sipros:latest'
    label 'sipros_large'    

    input:
    tuple path(config_file), path(ft_files)

    output:
    path 'sip/'

    script:
    """
    mkdir sip
    conda run -n sipros_env SiprosV4OMP -f *.FT2 -c $config_file -o sip/
    """
}

process sipros_PSM_filter {
    container 'stavisvols/psp_sipros:latest'
    label 'sipros_small'

    input:
    tuple path(config_file), path(sipfiles)

    output:
    tuple path(config_file), path('sip/')

    script:
    """
    conda run -n sipros_env python /opt/conda/envs/sipros_env/V4Scripts/sipros_peptides_filtering.py -c $config_file -w sip/
    """
}

process sipros_protein_filter {
    container 'stavisvols/psp_sipros:latest'
    label 'sipros_small'

    input:
    tuple path(config_file), path(sipfiles)

    output:
    tuple path(config_file), path('sip/')

    script:
    """
    conda run -n sipros_env python /opt/conda/envs/sipros_env/V4Scripts/sipros_peptides_assembling.py -c $config_file -w sip/
    """
}

process sipros_abundance_cluster {
    container 'stavisvols/psp_sipros:latest'
    label 'sipros_small'

    input:
    tuple path(config_file), path(sipfiles)

    output:
    path 'sip/'

    script:
    """
    conda run -n sipros_env python /opt/conda/envs/sipros_env/V4Scripts/ClusterSip.py -c $config_file -w sip/
    """
}

process sipros_protein_FDR {
    container 'stavisvols/psp_sipros:latest'
    label 'sipros_small'

    input:
    tuple path(sipfiles), val(row)

    output:
    tuple path('sip/'), val(row)

    script:
    """
    conda run -n sipros_env Rscript /opt/conda/envs/sipros_env/V4Scripts/refineProteinFDR.R -pro sip/*.pro.txt -psm sip/*.psm.txt -fdr 0.01 -o sip/$row.sample_ID
    """
}

process sipros_SIP_abundance {
    container 'stavisvols/psp_sipros:latest'
    publishDir path: "${params.results_dir}/${row.sample_ID}", mode: 'copy', pattern: "sip/*.*[!sip]"
    label 'sipros_med'

    input:
    tuple path(sipfiles), val(row)

    output:
    tuple path('sip/'), val(row)

    script:
    """
    conda run -n sipros_env Rscript /opt/conda/envs/sipros_env/V4Scripts/getLabelPCTinEachFT.R -pro sip/*.proRefineFDR.txt -psm sip/*.psm.txt -thr 3 -o sip/$row.sample_ID
    """
}

workflow sipros {
    take:
    row

    main:
    //set up per-file data as value channels
    row_channel = channel.value(row)
    config_file = channel.value(file(row.sipros_config))
    rawfile = channel.value(file(row.raw_file))
    ft_files = sipros_convert_raw_file(raw_file)
        | collect

    //run searches at each % RIA step
    config_files = sipros_config_generator(row_channel, config_file)
        | flatten
    search_results = sipros_search(config_files, ft_files)
        | collect
    
    //do post-processing
    processed_results = sipros_PSM_filter(config_file, search_results)
        | sipros_protein_filter
        | sipros_abundance_cluster
        | combine(row_channel)
        | sipros_protein_FDR
        | sipros_SIP_abundance

    emit:
    processed_results.out
}

workflow {
    //make results directory
    file(params.results_dir).mkdir()
    file(params.design).copyTo(params.results_dir)

    //parse the design file
    design = Channel.of(file(params.design)).splitCsv(header : true, sep : '\t', strip : true)

    //do per-file processing
    ipm_step_1 = sipros(design)
        // | sipros_to_ipm
        // | ipm_parse_data
        // | collect
    
    //run IPM classifier step

    //run IPM fitting jobs

    //
}