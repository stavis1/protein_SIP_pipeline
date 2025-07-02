
process config_generator {
    container 'stavisvols/psp_sipros:latest'
    label 'small'
    // publishDir params.results_dir, mode: 'copy'

    input:
    tuple val(row), path(global_config_file)

    output:
    tuple val(row), path('cfg/*.cfg'), path('config_*.cfg')

    script:
    """
    pathhash=\$(basename \$(pwd))
    sed -e '/IsopacketModeler/,\$d' $global_config_file > config_\$pathhash.cfg
    mkdir cfg
    elm=$row.label_elm
    if [ -z \$elm ]; then
        elm=C
    fi
    conda run -n sipros_env configGenerator -i config_\$pathhash.cfg -o cfg/ -e \$elm
    """
}

process convert_raw_file {
    container 'stavisvols/psp_sipros:latest'
    label 'med'

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

process search {
    container 'stavisvols/psp_sipros:latest'
    label 'large'    

    input:
    tuple val(sample_ID), path(ft_files), path(label_config_file), path(global_config_file), path(fasta)

    output:
    tuple val(sample_ID), path('*.sip'), path(global_config_file), path(fasta)

    script:
    """
    export OMP_NUM_THREADS=4
    conda run -n sipros_env SiprosV4OMP -f ./*.FT2 -c $label_config_file -o ./
    """
}

process psm_filter {
    container 'stavisvols/psp_sipros:latest'
    label 'small'

    input:
    tuple val(sample_ID), path(config_file), path(sipfiles), path(fasta)

    output:
    tuple val(sample_ID), path(config_file), path('sip/'), path(fasta)

    script:
    """
    mkdir sip
    cd sip
    ln -s ../*.sip .
    cd ../
    conda run -n sipros_env python /opt/conda/envs/sipros_env/V4Scripts/sipros_peptides_filtering.py -c $config_file -w sip/
    """
}

process protein_filter {
    container 'stavisvols/psp_sipros:latest'
    label 'small'

    input:
    tuple val(sample_ID), path(config_file), path(sipfiles), path(fasta)

    output:
    tuple val(sample_ID), path(config_file), path('sip/'), path(fasta)

    script:
    """
    conda run -n sipros_env python /opt/conda/envs/sipros_env/V4Scripts/sipros_peptides_assembling.py -c $config_file -w sip/
    """
}

process abundance_cluster {
    container 'stavisvols/psp_sipros:latest'
    label 'small'

    input:
    tuple val(sample_ID), path(config_file), path(sipfiles), path(fasta)

    output:
    tuple val(sample_ID), path('sip/'), path(fasta)

    script:
    """
    conda run -n sipros_env python /opt/conda/envs/sipros_env/V4Scripts/ClusterSip.py -c $config_file -w sip/
    """
}

process protein_FDR {
    container 'stavisvols/psp_sipros:latest'
    label 'small'

    input:
    tuple path(sipfiles), path(fasta), val(row)

    output:
    tuple path('sip/'), path(fasta), val(row)

    script:
    """
    conda run -n sipros_env Rscript /opt/conda/envs/sipros_env/V4Scripts/refineProteinFDR.R -pro sip/*.pro.txt -psm sip/*.psm.txt -fdr 0.01 -o sip/$row.sample_ID
    """
}

process sip_abundance {
    container 'stavisvols/psp_sipros:latest'
    publishDir path: "${params.results_dir}/${row.sample_ID}", mode: 'copy', pattern: "sip/*.txt"
    label 'med'

    input:
    tuple path(sipfiles), path(fasta), val(row)

    output:
    tuple val(row), path('sip/*.txt') 

    script:
    """
    conda run -n sipros_env Rscript /opt/conda/envs/sipros_env/V4Scripts/getLabelPCTinEachFT.R -pro sip/*.proRefineFDR.txt -psm sip/*.psm.txt -thr 3 -o sip/$row.sample_ID
    """
}

workflow sipros {
    take:
    rows

    main:
    //config file processing
    indexed_rows = rows.map {r -> tuple(r.sample_ID, r)}

    config_files = rows.map {r -> tuple(r, file(r.config))}
        | config_generator
        | flatMap {r -> r[1].collect {f -> tuple(r[0], f, r[2])}}
        | filter {row, cfg, cfg_g -> 
            def match = cfg.getBaseName() =~ /(\d+)Pct/
            def pct = match[0][1].toInteger()
            def reduction_factor = row.sipros_reduce.toInteger()*1000
            pct % reduction_factor == 0 || pct == 1070 
        }
        | map {row, cfg, cfg_g -> tuple(row.sample_ID, cfg, cfg_g)}

    search_jobs = rows.map {r -> tuple(r, file(r.raw_file))}
        | convert_raw_file
        //database searches are parallelized across config files
        | cross(config_files)
        | map {FT, config -> tuple(FT[0], FT[1], config[1], config[2])}

    search_results = rows.map {r -> tuple(r.sample_ID, file(r.fasta))}
        | cross(search_jobs)
        | map {fasta, job -> tuple(job[0], job[1], job[2], job[3], fasta[1])}
        | search
        | groupTuple(size: 101, remainder: true)
        | map {key, sips, configs, fastas -> tuple(key, configs[0], sips, fastas[0])}
        //post-processing
        | psm_filter
        | protein_filter
        | abundance_cluster
        | cross(indexed_rows)
        | map {abund, row -> tuple(abund[1], abund[2], row[1])}
        | protein_FDR
        | sip_abundance

    emit:
    search_results
}