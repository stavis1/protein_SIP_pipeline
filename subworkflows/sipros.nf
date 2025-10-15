
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
    label 'small'

    input:
    tuple val(row), path(rawfile)

    output:
    tuple val(row.sample_ID), path('*.FT{1,2}')

    script:
    //figure out how you want to determine the number of allocated cores then pass that to Raxport with the -j flag
    """
    mono /software/Sipros4/bin/Raxport.exe -i ./ -o ./ -j 3
    """
}

process process_fasta {
    container 'stavisvols/psp_sipros_python:latest'
    label 'med'

    input:
    tuple val(sample_ID), path(fasta), path(config), val(row)

    output:
    tuple val(sample_ID), path(fasta), path(config)

    script:
    """
    python /software/Sipros4/EnsembleScripts/sipros_prepare_protein_database.py -i $fasta -o DECOY_$fasta -c $config
    mv DECOY_$fasta fasta
    """
}

process unlabeled_search {
    container 'stavisvols/psp_sipros_mono:latest'
    label 'large'    

    input:
    tuple val(sample_ID), path(ft_files), path(fasta), path(config)

    output:
    tuple val(sample_ID), path('unlabeled_results/'), path(fasta), path(config)

    script:
    """
    export OMP_NUM_THREADS=4
    /software/Sipros4/bin/SiprosEnsembleOMP -f ./*.FT2 -c $config -o unlabeled_results/
    """
}

process unlabeled_python_postprocessing {
    container 'stavisvols/psp_sipros_python:latest'
    label 'small'

    input:
    tuple val(sample_ID), path(unlabeled_results), path(fasta), path(config)

    output:
    tuple val(sample_ID), path('unlabeled_results/'), path(fasta), path(config)

    script:
    """
    cp -rL unlabeled_results tmp
    rm unlabeled_results
    mv tmp unlabeled_results

    python /software/Sipros4/EnsembleScripts/sipros_psm_tabulating.py -i unlabeled_results/ -o unlabeled_results/ -c $config
    python /software/Sipros4/EnsembleScripts/sipros_ensemble_filtering.py -i unlabeled_results/ -o unlabeled_results/ -c $config
    python /software/Sipros4/EnsembleScripts/sipros_peptides_assembling.py -w unlabeled_results/ -c $config
    """
}
process unlabeled_R_postprocessing {
    container 'stavisvols/psp_sipros_r:latest'
    label 'small'

    input:
    tuple val(sample_ID), path(unlabeled_results), path(fasta), path(config)

    output:
    tuple val(sample_ID), path('unlabeled_results/'), path(fasta), path(config)

    script:
    """
    cp -rL unlabeled_results tmp
    rm unlabeled_results
    mv tmp unlabeled_results

    Rscript /software/Sipros4/V4Scripts/refineProteinFDR.R -pro $config/*.pro.txt -psm $config/*.psm.txt -fdr 0.01 -o $config/coli
    Rscript /software/Sipros4/V4Scripts/getSpectraCountInEachFT.R -pro $config/*.proRefineFDR.txt -psm $config/*.psm.txt -o $config/coli
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
    container 'stavisvols/psp_sipros_python:latest'
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
    python /software/Sipros4/V4Scripts/sipros_peptides_filtering.py -c $config_file -w sip/
    """
}

process protein_filter {
    container 'stavisvols/psp_sipros_python:latest'
    label 'small'

    input:
    tuple val(sample_ID), path(config_file), path(sipfiles), path(fasta)

    output:
    tuple val(sample_ID), path(config_file), path('sip/'), path(fasta)

    script:
    """
    python /software/Sipros4/V4Scripts/sipros_peptides_assembling.py -c $config_file -w sip/
    """
}

process abundance_cluster {
    container 'stavisvols/psp_sipros_python:latest'
    label 'small'

    input:
    tuple val(sample_ID), path(config_file), path(sipfiles), path(fasta)

    output:
    tuple val(sample_ID), path('sip/'), path(fasta)

    script:
    """
    python /software/Sipros4/V4Scripts/ClusterSip.py -c $config_file -w sip/
    """
}

process protein_FDR {
    container 'stavisvols/psp_sipros_r:latest'
    label 'small'

    input:
    tuple path(sipfiles), path(fasta), val(row)

    output:
    tuple path('sip/'), path(fasta), val(row)

    script:
    """
    Rscript /software/Sipros4/V4Scripts/refineProteinFDR.R -pro sip/*.pro.txt -psm sip/*.psm.txt -fdr 0.01 -o sip/$row.sample_ID
    """
}

process sip_abundance {
    container 'stavisvols/psp_sipros_r:latest'
    publishDir path: "${params.results_dir}/${row.sample_ID}", mode: 'copy', pattern: "sip/*.txt"
    label 'med'

    input:
    tuple path(sipfiles), path(fasta), val(row)

    output:
    tuple val(row), path('sip/*.txt') 

    script:
    """
    Rscript /software/Sipros4/V4Scripts/getLabelPCTinEachFT.R -pro sip/*.proRefineFDR.txt -psm sip/*.psm.txt -thr 3 -o sip/$row.sample_ID
    """
}

workflow sipros {
    take:
    rows

    main:
    //convert raw files to FT{1,2}
    converted_files = rows.map {r -> tuple(r, file(r.raw_file))}
        | convert_raw_file
    
    //add decoys to fasta files for unlabeled search
    fastas_with_decoys = rows.map {r -> tuple(r.sample_ID, file(r.fasta), file(r.config), r)}
        | process_fasta
        | branch { r -> 
            labeled: r.label_elm != null
            unlabeled: true
        }
    
    //unlabeled search
    unlabeled_search_results = fastas_with_decoys.unlabeled.cross(convert_raw_file)
        | map {fasta, FT -> tuple(FT[0], FT[1], fasta[1], fasta[2])}
        | unlabeled_search
        | unlabeled_python_postprocessing
        | unlabeled_R_postprocessing

    //config file processing for labeled search
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

    //set up list of individual atom % enrichment level search jobs
    search_jobs = converted_files.cross(config_files)
        | map {FT, config -> tuple(FT[0], FT[1], config[1], config[2])}

    //this allows rows to be added back for final processing 
    indexed_rows = rows.map {r -> tuple(r.sample_ID, r)}

    //run sipros SIP search on 
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