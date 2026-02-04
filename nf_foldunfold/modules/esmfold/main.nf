process SPLIT_MULTI {

    input:
    path (multi_fasta)

    output:
    path "*.fasta"

    script:
    """
    awk '/^>/ {if (f) close(f); f = substr(\$0, 2) ".fasta"} {print >> f}' $multi_fasta
    """
}

process RUN_ESMFOLD {

    conda '/path/to/envs/esmfold'

    input:
    tuple path(fasta), path(fold), path(model_dir)

    output:
    path "esmfold_$fasta.baseName"

    script:
    """
    mkdir esmfold_$fasta.baseName
    python3 ${fold} -i ${fasta} -o esmfold_$fasta.baseName --chunk-size 128 --model-dir ${model_dir} > esmfold_$fasta.baseName/log.txt
    """


    }
