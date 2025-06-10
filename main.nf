params.results_dir = "$launchDir/results"
//required columns in design file : sample_ID, label_elm, raw_file, sipros_config

process sipros_config_generator {
    container 'stavisvols/psp_sipros:latest'
    label 'sipros_small'
    // publishDir params.results_dir, mode: 'copy'

    input:
    tuple val(row), path(global_config_file)

    output:
    tuple val(row.sample_ID), path('cfg/*.cfg'), path(global_config_file)

    script:
    """
    mkdir cfg
    conda run -n sipros_env configGenerator -i $global_config_file -o cfg/ -e $row.label_elm
    """
}

process sipros_convert_raw_file {
    container 'stavisvols/psp_sipros:latest'
    label 'sipros_med'

    input:
    tuple val(row), path(rawfile)

    output:
    tuple val(row.sample_ID), path('*.FT{1,2}')

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
    tuple val(sample_ID), path(ft_files), path(label_config_file), path(global_config_file), path(fasta)

    output:
    tuple val(sample_ID), path('*.sip'), path(global_config_file)

    script:
    """
    export OMP_NUM_THREADS=4
    conda run -n sipros_env SiprosV4OMP -f ./*.FT2 -c $label_config_file -o ./
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
    mkdir sip
    cd sip
    ln -s ../*.sip .
    cd ../
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
    rows

    main:
    fasta_files = rows.map {r -> tuple(r.sample_ID, file(r.fasta))}

    //config file processing
    config_files = rows.map {r -> tuple(r, file(r.sipros_config))}
        | sipros_config_generator
        | flatMap {r -> r[1].collect {f -> tuple(r[0], f, r[2])}}

    search_results = rows.map {r -> tuple(r, file(r.raw_file))}
        | sipros_convert_raw_file
        //database searches are parallelized across config files
        | cross(config_files)
        | map {FT, config -> tuple(FT[0], FT[1], config[1], config[2])}
        | cross(fasta_files)
        | map {data, fasta -> tuple(data[0], data[1], data[2], data[3], fasta[1])}
        | sipros_search
        | groupTuple(size: 100)
        | map {k,v -> tuple(v[0][0], v[1].collect {e -> e[1]})}
        //post-processing
        | sipros_PSM_filter
        | sipros_protein_filter
        | sipros_abundance_cluster
        | combine(rows)
        | sipros_protein_FDR
        | sipros_SIP_abundance

    emit:
    search_results
}

workflow {
    //make results directory
    file(params.results_dir).mkdir()
    file(params.design).copyTo(params.results_dir)

    Channel.of(file(params.design)).splitCsv(header : true, sep : '\t', strip : true)
        | sipros
        // | collect
        // | isopacketModeler
        
}