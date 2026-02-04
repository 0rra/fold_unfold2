#!/bin/bash

out_dir="part1_refoldunfold/1_sequence_selection/sequences/swissprot_seqs"
script_dir="part1_refoldunfold/1_sequence_selection/scripts"

for i in 10 16 20 30 40 50 60 70 80 90 100 120 140 160 180 200;
# spliting file into files of seqs, specified seq length
do python ${script_dir}/split_fasta.py ${out_dir}/uniprot_sprot.fasta ${i} > ${out_dir}/${i}_seqs.fasta

# updating fasta headers to only contain uniprot id (no special characters)
python ${script_dir}/fix_names.py ${out_dir}/${i}_seqs.fasta > ${out_dir}/fixed_${i}_seqs.fasta

# running iupred to get sequence disorder scores
# may need updating to user iupred2a path
python scripts/run_iupred2a.py ${out_dir}/fixed_${i}_seqs.fasta ${out_dir}/fixed_${i}_seqs.results

# filtering fasta with iupred results to get non-disordered sequences
python ${script_dir}/filteriupred2a.py ${out_dir}/fixed_${i}_seqs.results
grep -A 1 -f ${out_dir}/passed_${i}_seqs ${out_dir}/fixed_${i}_seqs.fasta | sed '/^--/d' > ${out_dir}/filtered_${i}_seqs.fasta

# filter out fragment sequences
python ${script_dir}/filter_fragments.py ${out_dir}/filtered_${i}_seqs.fasta > ${out_dir}/filtered2_${i}_seqs.fasta

# joining back up randomly selected seqs
cat ${out_dir}/filtered2_${i}_seqs.fasta |awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' |shuf |head -n 10 |awk '{printf("%s\n%s\n",$1,$2)}' > ${out_dir}/sp_${i}_seqs.fasta
cat ${out_dir}/sp_${i}_seqs.fasta >> ${out_dir}/sp_seqs.fasta;
done


