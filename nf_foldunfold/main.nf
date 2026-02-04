nextflow.enable.dsl=2

include { RUN_COLABFOLD_BATCH } from "./modules/colabfold"
include { RUN_ESMFOLD } from "./modules/esmfold"
include { RUN_ALPHAFOLD3 } from "./modules/alphafold3"
include { FASTA_TO_JSON } from "./modules/alphafold3"
include { SPLIT_MULTI } from "./modules/esmfold"

workflow {
    main:

    if (params.input != null) {
        input_file = file(params.input)
    } else {
        input_file = null
    }
    // multi fasta file for colabfold
    multi_fasta = Channel.fromPath( input_file, checkIfExists: true )
    // single fasta channel for esmfold
    ch_fasta = SPLIT_MULTI(multi_fasta).flatMap()

    msa_flag  = Channel.value(params.use_msa)
    
    Channel.from(params.predictmethods.split(','))
    .branch { method ->

    esmfold: method == "esmfold"
        return [
        params.methods."${method}".fold,
        params.methods."${method}".model_dir
        ]

    colabfold: method == "colabfold"
        return [
        params.use_msa
        ]

    alphafold3: method == "alphafold3"
        return [
        params.methods."${method}".run_af3,
        params.methods."${method}".model_dir,
        params.methods."${method}".db_dir
        ]
    }.set { member_params }

    esm_params = ch_fasta.combine(member_params.esmfold)
    RUN_ESMFOLD(esm_params)
    colabfold_params = multi_fasta.combine(member_params.colabfold)
    RUN_COLABFOLD_BATCH(colabfold_params)
    fasta_params = multi_fasta.combine(msa_flag)
    json_dir = FASTA_TO_JSON(fasta_params)
    alphafold3_params = json_dir.combine(member_params.alphafold3)
    RUN_ALPHAFOLD3(alphafold3_params)
}

