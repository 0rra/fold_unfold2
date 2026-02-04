process FASTA_TO_JSON {

    input:
    tuple path(fasta), val(use_msa)

    output:
    path "af3_json_$fasta.baseName"

    script:
    """
    python $projectDir/scripts/fasta2json.py ${fasta} af3_json_$fasta.baseName ${use_msa}
    """

}


process RUN_ALPHAFOLD3 {

    input:
    tuple path(json_dir), path(run_af3), path(models), path(dbs)

    output:
    path "$json_dir.baseName"

    script:
    """
    mkdir -p $json_dir.baseName
    python ${run_af3} --input_dir=${json_dir} --model_dir=${models} --db_dir=${dbs} --output_dir=$json_dir.baseName
    """
}
