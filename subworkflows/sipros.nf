
process config_generator {
    container 'stavisvols/psp_sipros_mono:latest'
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
    /software/Sipros4/bin/configGenerator -i config_\$pathhash.cfg -o cfg/ -e \$elm
    """
}

process convert_raw_file {
    container 'stavisvols/psp_sipros_mono:latest'
    label 'medium'

    input:
    tuple val(row), path(rawfile)

    output:
    tuple val(row.sample_ID), path('*.FT{1,2}')

    script:
    //figure out how you want to determine the number of allocated cores then pass that to Raxport with the -j flag
    """
    if [[ ( $rawfile == *.raw ) || ( $rawfile == *.RAW ) ]]
    then
        mono /software/Sipros4/bin/Raxport.exe -i ./ -o ./ -j 3
    else
    if [[ ( $rawfile == *.zip ) ]]
    then
        unzip $rawfile
    fi
    fi
    """
}

process search {
    container 'stavisvols/psp_sipros_mono:latest'
    label 'large'    

    input:
    tuple val(sample_ID), path(ft_files), path(label_config_file), path(global_config_file), path(fasta)

    output:
    tuple val(sample_ID), path('*.sip'), path(global_config_file), path(fasta)

    script:
    """
    export OMP_NUM_THREADS=4
    /software/Sipros4/bin/SiprosV4OMP -f ./*.FT2 -c $label_config_file -o ./
    """
}

process psm_filter {
    container 'stavisvols/psp_percolator:latest'
    label 'percolator'

    input:
    tuple val(sample_ID), path(config_file), path(sipfiles), path(fasta)

    output:
    tuple val(sample_ID), path(config_file), path('*.txt'), path(fasta)

    script:
    """
    python /software/percolator.py $config_file
    """
}

process protein_filter {
    container 'stavisvols/psp_sipros_python:latest'
    label 'small'

    input:
    tuple val(sample_ID), path(config_file), path(psmfiles), path(fasta)
    
    output:
    tuple val(sample_ID), path(psmfiles), path('*.txt'), path(fasta)

    script:
    """
    python /software/Sipros4/V4Scripts/sipros_peptides_assembling.py -c $config_file -w ./
    python /software/Sipros4/V4Scripts/ClusterSip.py -c $config_file -w ./
    """
}

process sip_abundance {
    container 'stavisvols/psp_sipros_r:latest'
    publishDir path: "${params.results_dir}/${row.sample_ID}", mode: 'copy', pattern: "*.txt"
    label 'medium'

    input:
    tuple path(psmfiles), path(txtfiles), path(fasta), val(row)

    output:
    tuple val(row), path(psmfiles), path(txtfiles), path('*.txt') 

    script:
    """
    Rscript /software/Sipros4/V4Scripts/refineProteinFDR.R -pro *.pro.txt -psm *.psm.txt -fdr 0.01 -o $row.sample_ID
    Rscript /software/Sipros4/V4Scripts/getLabelPCTinEachFT.R -pro *.proRefineFDR.txt -psm *.psm.txt -thr 3 -o $row.sample_ID
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
        | cross(indexed_rows)
        | map {abund, row -> tuple(abund[1], abund[2], abund[3], row[1])}
        | sip_abundance

    emit:
    search_results
}