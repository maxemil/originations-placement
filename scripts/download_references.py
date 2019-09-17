import requests
import shutil
import os
import sys
import glob



gains = set()
for f in glob.glob('{}/*'.format(sys.argv[1])): # sys.argv[1] is a directory that contains tables from the ALE analysis (tables-gains-0.3)
    for line in open(f):
        if not line.startswith('cluster'):
            line = line.strip().split('\t')
            gains.add(line[0])

for g in gains:
    if not os.path.isfile('queries_split_all/{}.fasta'.format(g)):
        shutil.copyfile('../21_clusters_to_trees/clusters_eurNOG/faa/{}.fasta'.format(g), 'queries_split_all/{}.fasta'.format(g))

from collections import defaultdict

eur2NOG = defaultdict(list)
for line in open("../20_clusters_eggnog/NOG2EurNOG.tab"):
    line = line.strip().split('\t')
    for eur in line[1].split(';'):
        eur2NOG[eur].append(line[0])

cog2eur = defaultdict(list)
for f in glob.glob('queries_split_all/*'):
    f = f.replace('queries_split_all/', '').replace(".fasta", "")
    if len(f.split('_')) > 1:
        cog2eur[f.split('_')[1]].append(f)
    else:
        assert len(eur2NOG[f]) == 1
        cog2eur[eur2NOG[f][0]].append(f)

with open('COG2eurNOG_split.tsv', 'w') as out:
    for k,v in cog2eur.items():
        print("{}\t{}".format(k, ";".join(v)), file=out)


eggnog_refs = list(cog2eur.keys())

def download_from_eggnog(attribute, nogname, gzip, outdir):
  if not (os.path.isfile("{}/{}.{}.gz".format(outdir, nogname, attribute)) or os.path.isfile("{}/{}.{}".format(outdir, nogname, attribute))):
      url = "http://eggnogapi.embl.de/nog_data/file/{}/{}".format(attribute, nogname)
      r = requests.get(url, stream=True)
      if gzip:
         filename = "{}/{}.{}.gz".format(outdir, nogname, attribute)
      else:
         filename = "{}/{}.{}".format(outdir, nogname, attribute)
      with open(filename, 'wb') as f:
         shutil.copyfileobj(r.raw, f)

for eggnog in eggnog_refs:
    try:
        download_from_eggnog('tree', eggnog, True, 'references_trimal')
        download_from_eggnog('trimmed_alg', eggnog, True, 'references_trimal')
    except:
        print(eggnog)


import subprocess
import glob

for f in glob.glob("references_trimal/*.gz"):
    try:
        subprocess.check_call(["gunzip", f])
    except:
        subprocess.check_call(["mv", f, f.replace(".gz", "")])

from Bio import SeqIO
import os
import glob

taxids = [line.strip() for line in open("/local/two/Software/EggNOG_placement/Euryarchaeota.taxids")]

for ref in glob.glob("references_trimal/*.trimmed_alg"):
    if [rec.id.split('.')[0] in taxids for rec in SeqIO.parse(ref,'fasta')].count(False) < 4 \
        or len([rec for rec in SeqIO.parse(ref,'fasta')]) < 4 \
        or any([line.strip() == "Data not available" for line in open(ref)]):
       os.system("mv {} failed_references_trimal".format(ref))
       os.system("mv {} failed_references_trimal".format(ref.replace('trimmed_alg', 'tree')))
