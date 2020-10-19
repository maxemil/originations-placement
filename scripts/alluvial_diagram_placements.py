#!/usr/bin/env python3
"""
  19-04-04, Max Schön, <max-emil.schon@icm.uu.se>
  example usage:
"""
import os
import argparse
import pandas as pd
import glob
from ete3 import ncbi_taxonomy
from Bio import SeqIO
from collections import defaultdict
import subprocess

def get_phylum(taxon):
    # # ref_taxa = set()
    ref_taxa = set(['Proteobacteria', 'Euryarchaeota', 'Terrabacteria group',
                    'Acidobacteria', 'Spirochaetes', 'TACK group', 'PVC group',
                    'FCB group', 'Deferribacteres', 'Aquificae', 'DPANN group',
                    'Synergistetes', 'Eukaryota', 'Thermotogae'])
                    # 'Metazoa', 'Fungi', 'Viridiplantae', 'Alveolata',
                    # 'Rhodophyta', 'Fornicata', 'Choanoflagellida',
                    # 'Amoebozoa', 'Stramenopiles'])
    others = set(['Acidobacteria', 'Aquificae', 'Deferribacteres',
                  'DPANN group', 'PVC group', 'Synergistetes', 'Thermotogae'])
    try:
        lineage = ncbi.get_lineage(ncbi.get_name_translator([taxon])[taxon][0])
        for l in lineage:
            lname = ncbi.get_taxid_translator([l])[l]
            if lname in ref_taxa:
                if lname in others:
                    return 'others'
                return lname
            # if ncbi.get_rank([l])[l] == 'phylum':
            #     return ncbi.get_taxid_translator([l])[l]
        return taxon
    except:
        return taxon

def print_csv_plot_alluvial(node, df):
    df = df.drop("node", axis=1)
    df.to_csv('{}_freqs.csv'.format(node), sep='\t', index=False)
    with open("alluvial_plot.R", 'w') as out:
        print("library(alluvial)", file=out)
        print("placement = read.table('{}_freqs.csv', sep='\t', header=T)".format(node), file=out)
        print("pdf('alluvial_{}.pdf')".format(node), file=out)
        #hide = placement$tax == '',
        print("    alluvial(placement[,1:3], freq=placement$freq, cex=0.5, col=ifelse(placement$event == 'gain', 'orange', 'grey' ))", file=out)
        print("dev.off()", file=out)
    subprocess.call("Rscript alluvial_plot.R".split())

# ordered_tax = c('de-novo genes', 'others', 'cellular organisms', 'Eukaryota',
#                 'Archaea', 'Euryarchaeota', 'TACK group', 'DPANN group',
#                 'Bacteria', 'Proteobacteria', 'Terrabacteria group',
#                 'Acidobacteria', 'Spirochaetes', 'PVC group',
#                 'FCB group', 'Deferribacteres', 'Aquificae',
#                 'Synergistetes', 'Thermotogae', '',
#                 'Unknown', 'Metabolism', 'Cellular Processes', 'Information',
#                 'gain', 'loss')
#
# placement = read.table('81_freqs.csv', sep='\t', header=T)
# placement = as.data.frame(placement)
# # df.lode = to_lodes_form(placement, axes = 1:2, id=id)
# # df.lode = mutate(df.lode, level = factor(level, levels=ordered_tax))
# df.lode = mutate(placement, subject = seq(1, n())) %>% gather(x, level, -freq, -subject) %>% mutate(level = factor(level, levels=ordered_tax))
#
# ggplot(df.lode,
#        aes(x = x,
#            stratum = level,
#            alluvium = subject,
#            y = freq)) +
#   geom_alluvium(width = 1/12) +
#   geom_stratum(width = 1/12, fill = "white", color = "grey") +
#   geom_label(stat = "stratum", label.strata = TRUE, label.size = NA, fill = NA) +
#   ggtitle("taxonomic origin at node 81")
#   scale_x_discrete(limits = c("Category", "Taxonomy"), expand = c(.05, .05)) +
#
# ggplot(as.data.frame(placement), aes(y = freq, axis1 = tax, axis2 = cat, fill = event, group = factor(tax, levels=ordered_tax))) +
#   geom_alluvium(width = 1/12) +
#   geom_stratum(width = 1/12, fill = "white", color = "grey") +
#   scale_x_discrete(limits = c("Taxonomy", "Category"), expand = c(.05, .05)) +
#   geom_label(stat = "stratum", label.strata = TRUE, label.size = NA, fill = NA) +
#   ggtitle("taxonomic origin at node 81")

parser = argparse.ArgumentParser()

parser.add_argument("placements", type=str,nargs='+',
                    help="gappa csv output files")
parser.add_argument("-a", "--annotations", type=str,
                    help="EggNOG OG annotation file")

args = parser.parse_args()
ncbi = ncbi_taxonomy.NCBITaxa()


# annot = {line.split('\t')[0]:";".join(eval(line.split('\t')[4])) for line in open(args.annotations)}
annot = {line.split('\t')[0]:";".join(eval(line.split('\t')[4])) for line in open("eurNOG_annotations.tsv")}

tax = pd.DataFrame()
# for placement in args.placements:
for placement in glob.glob("placements/*.csv"):
    base = os.path.basename(placement).split('.')[0]
    cog = os.path.basename(placement).split('.')[1]
    df = pd.read_csv(placement, sep='\t')
    taxon = df['taxopath'][df['fract'].idxmax()].split(';')[-1]
    tax = tax.append(pd.DataFrame({"COG":cog, "tax":taxon}, index=[base]))

taxids = [line.strip() for line in open("/local/two/Software/EggNOG_placement/Euryarchaeota.taxids")]
eur2NOG = {}
for line in open('COG2eurNOG_split.tsv'):
    line = line.strip().split()
    nog = line[0]
    for eur in line[1].split(';'):
        eur2NOG[eur] = nog

for f in glob.glob("queries_split_all/*"):
    eur = f.replace('queries_split_all/', '').replace('.fasta', '')
    if not os.path.exists('placements/{}.{}.csv'.format(eur, eur2NOG[eur])):
        try:
            if [rec.id.split('.')[0] in taxids for rec in SeqIO.parse("failed_references/{}.raw_alg".format(eur2NOG[eur]),'fasta')].count(False) < 4:
                tax = tax.append(pd.DataFrame({"COG":eur2NOG[eur], "tax":"de-novo genes"}, index=[eur]))
        except:
            print("something else is wrong with reference {} of {}".format(eur2NOG[eur], eur))

clst2node = {}
node2tab = {}
for t in glob.glob('tables-gains-0.3/*.tab'):
    node = t.strip('.tab').split('to')[1]
    df = pd.read_csv(t, sep='\t')
    node2tab[node] = df
    for c in df['cluster']:
        if not c in clst2node.keys():
            clst2node[c] = node
        else:
            freq_curr = node2tab[node][node2tab[node]['cluster'] == c]['acquisition_freq'].item()
            freq_prev = node2tab[clst2node[c]][node2tab[clst2node[c]]['cluster'] == c]['acquisition_freq'].item()
            if float(freq_curr) > float(freq_prev):
                clst2node[c] = node
            else:
                pass
                # print("{} occurs more than once, currenct is {} at {}, prev was {} at {}".format(c, freq_curr, node, freq_prev, clst2node[c]))

loss = pd.DataFrame()
for t in glob.glob('tables-loss-0.3/*.tab'):
    node = t.strip('.tab').split('to')[1]
    n2cat = defaultdict(int)
    for line in open(t):
        if not line.startswith('cluster'):
            cat = line.split('\t')[2]
            n2cat[cat] += 1
    for k,v in n2cat.items():
        loss = loss.append(pd.DataFrame({"cat":k, "freq":v, "node":node}, index=[0]))


annot = pd.DataFrame.from_dict(annot, orient="index", columns=["cat"])
node = pd.DataFrame.from_dict(clst2node, orient="index", columns=["node"])

cat_groups = {"Metabolism":['C','G','E','F','H','I','P','Q', 'C;E', 'E;P'],
              "Cellular Processes":['D','Y','V','T','M','N','Z','W','U','O'],
              "Information":['J','A','K','L','B'],
              "Unknown":['S']}
catmap = defaultdict(lambda: "Unknown")
for k, value in cat_groups.items():
    for v in value:
        catmap[v] = k
annot['cat'] = annot['cat'].apply(lambda x: catmap[x])
loss['cat'] = loss['cat'].apply(lambda x: catmap[x])

df = tax.merge(annot, left_index=True, right_index=True, how='left').merge(node, left_index=True, right_index=True, how='left').fillna('Unknown')
df['tax'] = df['tax'].apply(get_phylum)
freq = df.groupby(['node', 'tax','cat']).size().reset_index()
freq.columns = ['node', 'tax', 'cat', 'freq']
freq["event"] = "gain"
loss["event"] = "loss"
freq = pd.concat([freq, loss], sort=True)
freq = freq[['tax', 'cat', 'event', 'freq', 'node']]
freq.reset_index(drop=True, inplace=True)
freq.to_csv('placement_freq.csv', sep='\t', index=False)
for g in freq.groupby(['node']):
    alluvial = g[1]
    if (g[1]['event'] == 'loss').sum() > (g[1]['event'] == 'gain').sum():
        alluvial = g[1][g[1]['event'] == 'gain']
    print_csv_plot_alluvial(g[0], alluvial)
