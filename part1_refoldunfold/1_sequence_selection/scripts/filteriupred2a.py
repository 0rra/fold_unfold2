import sys

iupred_res = sys.argv[1]

with open(iupred_res) as pred:
    seqpreds = {}
    seqid = ""
    for line in pred:
        # fix output format issue
        if line.startswith("# Results for"):
            id_line = line.strip()
            seqid = id_line.split("Results for")[-1]
            seqid = seqid.strip()
            seqpreds[seqid] = {}

        if line[0].isdigit():
            # print(line)
            info = line.split("\t")
            position = info[0]
            res = info[1]
            score = info[2]
            seqpreds[seqid][position] = float(score)

to_files = {}

# calc average
for i in seqpreds:
    sum = 0
    for res in seqpreds[i]:
        sum += seqpreds[i][res]
    length = len(seqpreds[i])
    average = sum / length
    seqpreds[i]["average"] = average
    if average <= 0.5:
        try:
            to_files[f"{length}"].append(i)
        except KeyError:
            to_files[f"{length}"] = [i]
        # open file

for l in to_files:
    with open(f"passed_{l}_seqs", "w") as resfile:
        resfile.write("\n".join(to_files[l]))
