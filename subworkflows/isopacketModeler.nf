
process sipros_psm_converter {
    container 'stavisvols/psp_isopacketmodeler'
    label 'ipm_small'

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
    container 'stavisvols/psp_isopacketmodeler'
    label 'ipm_3core'

    input:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer)
    
    output:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path('design.tsv'), path('*step1_*.dill')

    script:
    """
    echo -e 'file\\tlabel' > design.tsv
    if [ -z ${label_elm} ]; then
        echo -e '${mzml}\\t' >> design.tsv
    else
        echo -e '${mzml}\\t${label_elm}[${label_integer}]' >> design.tsv
    fi
    conda run -n isotope_env -m isopacketModeler \\
        --working_directory ./ \\
        --output_directory ./ \\
        --design_file design.tsv \\
        --mzml_dir ./ \\    
        --psms ${psms} \\
        --psm_headers seq,file,ScanNumber,ParentCharge,proteins \\
        --AA_formulae ${amino_acids} \\
        --cores 3 \\
        --data_generating_processes BetabinomQuiescentMix \\
        --data_generating_processes Betabinom \\
        --data_generating_processes BinomQuiescentMix \\
        --data_generating_processes Binom \\
        --do_PSM_classification \\
        --stopping_point 1
    """
}

process classifier {
    container 'stavisvols/psp_isopacketmodeler'
    label 'ipm_small'

    input:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(design_file), path(checkpoints)

    output:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(design_file), path('*step2_*.dill')

    script:
    """
    conda run -n isotope_env -m isopacketModeler -o ${options} \\
        --working_directory ./ \\
        --output_directory ./ \\
        --design_file design.tsv \\
        --mzml_dir ./ \\
        --psms ${psms} \\
        --psm_headers seq,file,ScanNumber,ParentCharge,proteins \\
        --AA_formulae ${amino_acids} \\
        --cores 3 \\
        --data_generating_processes BetabinomQuiescentMix \\
        --data_generating_processes Betabinom \\
        --data_generating_processes BinomQuiescentMix \\
        --data_generating_processes Binom \\
        --do_PSM_classification \\
        --checkpoint_files *step1_*.dill \\
        --stopping_point 2
    """
}

process scatter_peptides {
    container 'stavisvols/psp_isopacketmodeler'
    label 'ipm_small'

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
    container 'stavisvols/psp_isopacketmodeler'
    label 'ipm_small'

    input:
    tuple path(psms), path(mzml), path(amino_acids), val(label_elm), val(label_integer), path(design_file), path(checkpoints)

    output:
    tuple path('*peptides.dill'), path('*peptides.tsv')

    script:
    """
    conda run -n isotope_env -m isopacketModeler -o ${options} \\
        --working_directory ./ \\
        --output_directory ./ \\
        --design_file design.tsv \\
        --mzml_dir ./ \\
        --psms ${psms} \\
        --psm_headers seq,file,ScanNumber,ParentCharge,proteins \\
        --AA_formulae ${amino_acids} \\
        --cores 3 \\
        --data_generating_processes BetabinomQuiescentMix \\
        --data_generating_processes Betabinom \\
        --data_generating_processes BinomQuiescentMix \\
        --data_generating_processes Binom \\
        --do_PSM_classification \\
        --checkpoint_files subset_*.dill
    
    mv peptides.dill \$\$_peptides.dill
    mv peptides.tsv \$\$_peptides.tsv
    """
}

process merge_results {
    container 'stavisvols/psp_isopacketmodeler'
    label 'ipm_small'

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
        | collect
        | classifier
        | scatter_peptides
        | flatMap {psms, mzml, aas, labelE, labelI, design, dills -> dills.collect {dill -> tuple(psms, mzml, aas, labelE, labelI, design, dill)}}
        | model_fitting
        | collect
        | merge_results

    emit:
    peptides
}
