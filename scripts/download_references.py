import requests
import shutil
import os
import sys


eggnog_refs = [line.strip() for line in open(sys.argv[1])]

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
        download_from_eggnog('tree', eggnog, True, 'references')
        download_from_eggnog('raw_alg', eggnog, True, 'references')
    except:
        print(eggnog)


import subprocess
import glob

for f in glob.glob("references/*.gz"):
    try:
        subprocess.check_call(["gunzip", f])
    except:
        subprocess.check_call(["mv", f, f.replace(".gz", "")])

from Bio import SeqIO
import os
import glob

taxids = [line.strip() for line in open(sys.argv[2])]

for ref in glob.glob("references/*.raw_alg"):
    if [rec.id.split('.')[0] in taxids for rec in SeqIO.parse(ref,'fasta')].count(False) < 4 \
        or len([rec for rec in SeqIO.parse(ref,'fasta')]) < 4 \
        or any([line.strip() == "Data not available" for line in open(ref)]):
       os.system("mv {} failed_references".format(ref))
       os.system("mv {} failed_references".format(ref.replace('raw_alg', 'tree')))
