
process sipros_psm_converter {
    container 'stavisvols/psp_isopacketmodeler:latest'
    label 'small'

    input:
    tuple val(row), path(sipros_results)

    output:
    tuple val(row), path('*.psms')

    script:
    """
    python /scripts/sipros2IPM.py ${row.sample_ID}
    """
}

process convert_raw_file {
    container 'stavisvols/psp_sipros_mono:latest'
    label 'small'
    stageInMode 'link'

    input:
    tuple path(psms), path(raw_file), val(label_elm), val(label_integer), path(config)

    output:
    tuple path(psms), path('*.mzML'), val(label_elm), val(label_integer), path(config)

    script:
    """
    if [[ ( $raw_file == *.raw ) || ( $raw_file == *.RAW ) ]]
    then
        (timeout 10m mono /software/ThermoRawFileParser.exe -i $raw_file -o ./ -f 2; exit 0)
        ls *.mzML
    fi
    """
}

process parse_mzml_files {
    container 'stavisvols/psp_isopacketmodeler:latest'
    label 'med'

    input:
    tuple path(psms), path(mzml), val(label_elm), val(label_integer), path(config)
    
    output:
    tuple path(psms), path(mzml), path('AA_formulae.tsv'), val(label_elm), val(label_integer), path(config), path('design.tsv'), path('*step1_*.dill')

    script:
    """
    #parse config file
    psm_line=\$(echo "['"$psms"']")
    sed -n -e '/IsopacketModeler/,\$p' $config | \\
        sed 's|@CORES|$task.cpus|g' | \\
        sed "s|@PSMS|\$psm_line|g" | \\
        sed 's|@CHECKPOINT|\\[\\]|g' | \\
        sed 's|@STOP|1|g' > options.toml

    #make design file
    mzml=$mzml
    filename="\${mzml%.*}"
    echo -e 'file\\tlabel' > design.tsv 
    if [ -z ${label_elm} ]; then
        echo -e \$filename'\\t' >> design.tsv
    else
        echo -e \$filename'\\t${label_elm}[${label_integer}]' >> design.tsv
    fi

    #get amino acid formula file
    python /scripts/formula_parser.py $config

    #parse mzML file
    python -m isopacketModeler file --options options.toml
    pathhash=\$(basename \$(pwd))
    dillfile=\$(ls *step1_*.dill)
    mv \$dillfile \$(echo \$dillfile | sed "s|\\(.*_step1_\\)[[:digit:]]*.dill|\\1\${pathhash}.dill|g")
    """
}

process classifier {
    container 'stavisvols/psp_isopacketmodeler:latest'
    label 'huge'

    input:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(config), path(design_file), path(checkpoints)

    output:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(config), path(design_file), path('*step2_*.dill')

    script:
    """
    #parse config file
    psm_line=\$(echo "['"$psms"']")
    sed -n -e '/IsopacketModeler/,\$p' $config | \\
        sed 's|@CORES|$task.cpus|g' | \\
        sed "s|@PSMS|\$psm_line|g" | \\
        sed 's|@CHECKPOINT|\\["*step1_*.dill"\\]|g' | \\
        sed 's|@STOP|2|g' > options.toml
    python -m isopacketModeler file --options options.toml
    """
}

process scatter_peptides {
    container 'stavisvols/psp_isopacketmodeler:latest'
    label 'small'

    input:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(config), path(design_file), path(checkpoints)

    output:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(config), path(design_file), path('subset_*.dill')

    script:
    """
    python /scripts/split_peptides.py ${checkpoints} 50
    """
}

process model_fitting {
    container 'stavisvols/psp_isopacketmodeler:latest'
    label 'large'

    input:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(config), path(design_file), path(checkpoints)

    output:
    tuple path('*peptides.dill'), path('*peptides.tsv')

    script:
    """
    #parse config file
    psm_line=\$(echo "['"$psms"']")
    sed -n -e '/IsopacketModeler/,\$p' $config | \\
        sed 's|@CORES|$task.cpus|g' | \\
        sed "s|@PSMS|\$psm_line|g" | \\
        sed 's|@CHECKPOINT|\\["subset_*.dill"\\]|g' | \\
        sed 's|@STOP|false|g' > options.toml
    python -m isopacketModeler file --options options.toml
    
    pathhash=\$(basename \$(pwd))
    mv peptides.dill \${pathhash}_peptides.dill
    mv peptides.tsv \${pathhash}_peptides.tsv
    """
}

process merge_results {
    container 'stavisvols/psp_isopacketmodeler:latest'
    publishDir path: params.results_dir, mode: 'copy', pattern: "peptides.*"
    label 'small'

    input:
    tuple path(dill), path(tsv) 

    output:
    tuple path('peptides.dill'), path('peptides.tsv')

    script:
    """
    python /scripts/merge_peptides.py
    """
}


workflow sipros_to_ipm {
    take:
    samples //tuple(rows, sipros_results)

    main:
    psms = sipros_psm_converter(samples)

    emit:
    psms
}

workflow isopacketModeler {
    take:
    samples //tuple(rows, PSMs)

    main:
    peptides = samples.map {row, psms -> tuple(psms, file(row.raw_file), row.label_elm, row.label_integer, file(row.config))}
        | convert_raw_file
        | parse_mzml_files
        | collect(flat:false)
       	| map {data -> tuple(data[0][0..6] + [data.collect {f -> f[7]}])}
        | classifier
        | scatter_peptides
        | flatMap {psms, mzml, aas, labelE, labelI, config, design, dills -> dills.collect {dill -> tuple(psms, mzml, aas, labelE, labelI, config, design, dill)}}
        | model_fitting
        | collect
        | map {data -> [data.findAll {it.getExtension() == 'dill'}] + [data.findAll {it.getExtension() == 'tsv'}]}
        | merge_results

    emit:
    peptides
}
