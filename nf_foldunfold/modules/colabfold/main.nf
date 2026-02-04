process RUN_COLABFOLD_BATCH {

    conda '/path/to/env/localcolabfold/colabfold-conda'

    input:
    tuple path(fasta), val(use_msa)

    output:
    path "colabfold_$fasta.baseName"
    
    script:
    def msa_flag = use_msa ? "" : "--msa-mode single_sequence"

    """
    colabfold_batch --model-type auto --amber --num-relax 1 ${msa_flag} ${fasta} colabfold_${fasta.baseName}
    """

    }
