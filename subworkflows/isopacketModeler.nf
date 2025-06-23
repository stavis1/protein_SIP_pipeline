
process sipros_psm_converter {
    container 'stavisvols/psp_isopacketmodeler:latest'
    label 'small'

    input:
    tuple val(row), path(sipros_results)

    output:
    tuple val(row), path('*.psms')

    script:
    """
    conda run -n isotope_env python /scripts/sipros2IPM.py ${row.sample_ID}
    """
}


process parse_mzml_files {
    container 'stavisvols/psp_isopacketmodeler:latest'
    label 'med'

    input:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer)
    
    output:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path('design.tsv'), path('*step1_*.dill')

    script:
    """
    filename="\${${mzml}%.*}"
    echo -e 'file\\tlabel' > design.tsv 
    if [ -z ${label_elm} ]; then
        echo -e '\$filename\\t' >> design.tsv
    else
        echo -e '\$filename\\t${label_elm}[${label_integer}]' >> design.tsv
    fi
    cmd='''
    --working_directory ./
    --output_directory ./
    --design_file design.tsv
    --mzml_dir ./
    --psms ${psms}
    --psm_headers seq,file,ScanNumber,ParentCharge,proteins
    --aa_formulae ${amino_acids}
    --cores 3
    --data_generating_processes BetabinomQuiescentMix
    --data_generating_processes Betabinom
    --data_generating_processes BinomQuiescentMix
    --data_generating_processes Binom
    --do_psm_classification
    --stopping_point 1
    --overwrite
    '''
    conda run -n isotope_env python -m isopacketModeler cmd \$(echo \$cmd | tr -d '\\n')
    pathhash=\$(basename \$(pwd))
    dillfile=\$(ls *step1_*.dill)
    mv \$dillfile \$(echo \$dillfile | sed "s|\\(.*_step1_\\)[[:digit:]]*.dill|\\1\${pathhash}.dill|g")
    """
}

process classifier {
    container 'stavisvols/psp_isopacketmodeler:latest'
    label 'huge'

    input:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(design_file), path(checkpoints)

    output:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(design_file), path('*step2_*.dill')

    script:
    """
    cmd='''
    --working_directory ./
    --output_directory ./
    --design_file design.tsv
    --mzml_dir ./
    --psms ${psms}
    --psm_headers seq,file,ScanNumber,ParentCharge,proteins
    --aa_formulae ${amino_acids}
    --cores 3
    --data_generating_processes BetabinomQuiescentMix
    --data_generating_processes Betabinom
    --data_generating_processes BinomQuiescentMix
    --data_generating_processes Binom
    --do_psm_classification
    --stopping_point 2
    --overwrite
    '''
    conda run -n isotope_env python -m isopacketModeler cmd \$(echo \$cmd | tr -d '\\n') --checkpoint_files '*step1_*.dill'
    """
}

process scatter_peptides {
    container 'stavisvols/psp_isopacketmodeler:latest'
    label 'small'

    input:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(design_file), path(checkpoints)

    output:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(design_file), path('subset_*.dill')

    script:
    """
    conda run -n isotope_env python /scripts/split_peptides.py ${checkpoints} 50
    """
}

process model_fitting {
    container 'stavisvols/psp_isopacketmodeler:latest'
    label 'large'

    input:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(design_file), path(checkpoints)

    output:
    tuple path('*peptides.dill'), path('*peptides.tsv')

    script:
    """
    cmd='''
    --working_directory ./
    --output_directory ./
    --design_file design.tsv
    --mzml_dir ./
    --psms ${psms}
    --psm_headers seq,file,ScanNumber,ParentCharge,proteins
    --aa_formulae ${amino_acids}
    --cores 3
    --data_generating_processes BetabinomQuiescentMix
    --data_generating_processes Betabinom
    --data_generating_processes BinomQuiescentMix
    --data_generating_processes Binom
    --do_psm_classification
    --checkpoint_files subset_*.dill
    --overwrite
    '''
    conda run -n isotope_env python -m isopacketModeler cmd \$(echo \$cmd | tr -d '\\n')
    
    pathhash=\$(basename \$(pwd))
    mv peptides.dill \${pathhash}_peptides.dill
    mv peptides.tsv \${pathhash}_peptides.tsv
    """
}

process merge_results {
    container 'stavisvols/psp_isopacketmodeler:latest'
    label 'small'

    input:
    tuple path(dill), path(tsv) 

    output:
    tuple path('peptides.dill'), path('peptides.tsv')

    script:
    """
    conda run -n isotope_env python /scripts/merge_peptides.py
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
    peptides = samples.map {row, psms -> tuple(psms, file(row.mzml), file(row.amino_acids), row.label_elm, row.label_integer)}
        | parse_mzml_files
        | collect(flat:false)
       	| map {data -> tuple(data[0][0..5] + [data.collect {f -> f[6]}])}
        | classifier
        | scatter_peptides
        | flatMap {psms, mzml, aas, labelE, labelI, design, dills -> dills.collect {dill -> tuple(psms, mzml, aas, labelE, labelI, design, dill)}}
        | model_fitting
        | collect
        | merge_results

    emit:
    peptides
}
