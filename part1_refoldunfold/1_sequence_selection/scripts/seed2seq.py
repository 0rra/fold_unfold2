# Opens seed file and gets first 3 sequences for each family
with open("part1_refoldunfold/1_sequence_selection/sequences/antifam_seed_seqs/AntiFam.seed") as seed:
    seqs = {}
    family_id = ""
    found_gf_sq = False
    seq_lim = 3

    for line in seed:
        if line.startswith("#=GF AC"):
            family_id = line.replace("#=GF AC", "").strip()
            seqs[family_id] = []
            seq_lim = 3 

        if line.startswith("#=GF SQ"):
            seed_num = line.replace("#=GF SQ", "").strip()
            seed_num = int(seed_num)
            if seed_num < 3:
                seq_lim = seed_num
                print(f"Family {family_id} has only {seq_lim} sequences")
            found_gf_sq = True
            continue

        if found_gf_sq and family_id and len(seqs[family_id]) < seq_lim:
            if line.strip() and not line.startswith("#"):
                sequence = line.strip()
                sequence = sequence.split()[-1]
                sequence = sequence.replace("-", "").replace(".", "")

                if sequence:
                    seqs[family_id].append(sequence.upper())

                    if len(seqs[family_id]) == seq_lim:
                        found_gf_sq = False


n = 0
with open(
    "part1_refoldunfold/1_sequence_selection/sequences/antifam_seed_seqs/seed_seqs.fasta2", "w"
) as ff:
    for family_id in seqs:
        for count, sequence in enumerate(seqs[family_id], start=1):
            n += 1
            ff.write(f">{family_id}_seq_{count}\n")
            ff.write(f"{sequence}\n")

print(f"Written {n} sequences from {len(seqs)} families")
