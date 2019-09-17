from collections import defaultdict

for line in open("../20_clusters_eggnog/NOG2EurNOG.tab"):
    line = line.strip().split('\t')
    for eur in line[1].split(';'):
        eur2NOG[eur].append(line[0])

cog2eur = defaultdict(list)
for f in glob.glob('queries_split/*'):
    f = f.replace('queries_split/', '').replace(".fasta", "")
    if len(f.split('_')) > 1:
        cog2eur[f.split('_')[1]].append(f)
    else:
        assert len(eur2NOG[f]) == 1
        cog2eur[eur2NOG[f][0]].append(f)

with open('COG2eurNOG_split.tsv', 'w') as out:
    for k,v in cog2eur.items():
        print("{}\t{}".format(k, ";".join(v)), file=out)
