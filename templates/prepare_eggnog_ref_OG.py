#! /opt/conda/envs/placement/bin/python3

from Bio import SeqIO, Seq
from ete3 import ncbi_taxonomy, Tree
from itertools import combinations
ncbi = ncbi_taxonomy.NCBITaxa()

def commonprefix(taxon_paths):
    taxon_path_lists = [t.split('; ') for t in taxon_paths if t is not None]
    s1 = min(taxon_path_lists)
    s2 = max(taxon_path_lists)
    for i, c in enumerate(s1):
        if c != s2[i]:
            return s1[:i]
    return s1

def get_lineage(taxid):
    try:
        lineage = ncbi.get_lineage(int(taxid))
        lineage_string = ""
        for l in lineage:
            rank = ncbi.get_rank([l])[l]
            lineage_string = "{} {};".format(lineage_string, ncbi.get_taxid_translator([l])[l])
        return lineage_string
    except:
        return None

def clean_rec_id(recid):
    illegal_chars = ":,();[]' /+%\$@#&*<=>^!?\\{\\}\\"`"
    for c in illegal_chars:
        recid = recid.replace(c, '')
    return recid

def remove_ident_seqs_tax(node, seqdict, counter):
    nodes = node.children
    node_names = [n.name for n in nodes]
    #assert all([seqdict[n.name].seq == seqdict[p.name].seq for n,p in combinations(nodes, 2)])
    tax_strings = [get_lineage(n.name.split('.')[0]) for n in nodes]
    tax_path = commonprefix(tax_strings)
    taxon = tax_path[-1].replace(';', '')
    id = ncbi.get_name_translator([taxon])[taxon][0]
    node.name = "{}.{}_{}".format(id, taxon.replace(' ', '_'), counter)
    node.name = clean_rec_id(node.name)
    remove_ident_seqs_aln(node.name, node_names, seqdict)
    node.children = []

def remove_ident_seqs_aln(node_name, children, seqdict):
    seqdict[node_name] = seqdict[children[0]]
    seqdict[node_name].id = node_name
    seqdict[node_name].description = node_name
    seqdict[node_name].name = node_name
    for c in children:
        print(c)
        seqdict.pop(c)

def remove_polytomies(tree, seqdict):
    counter=0
    for node in tree.traverse(strategy='postorder'):
        if len(node.children) > 2 and not node.is_root():
            remove_ident_seqs_tax(node, seqdict, counter)
            counter += 1

def remove_taxids(tree, seqdict, taxid_file):
    taxids = [line.strip() for line in open(taxid_file)]
    keep_leaves = [l.name for l in tree.iter_leaves()]
    for l in tree.iter_leaves():
        if l.name.split('.')[0] in taxids:
            keep_leaves.remove(l.name)
            seqdict.pop(l.name)
    tree.prune(keep_leaves, preserve_branch_length=True)


def main(tree, aln, treeout, alnout, taxout, taxid_file, precompiled_taxa):
  tree = Tree(tree, format=0)
  seqdict = {rec.id: rec for rec in SeqIO.parse(aln, 'fasta')}
  remove_taxids(tree, seqdict, taxid_file)
  remove_polytomies(tree, seqdict)
  precompiled_taxa = {line.split('\\t')[0]:line.strip().split('\\t')[1] for line in open(precompiled_taxa)}
  with open(taxout, 'w') as out:
      for rec in seqdict.values():
          try:
              lineage = precompiled_taxa[rec.id.split('.')[0]]
          except:
              lineage = get_lineage(rec.id.split('.')[0])
          print("{}\\t{}".format(rec.id, lineage), file=out)
  with open(alnout, 'w') as out:
      SeqIO.write(seqdict.values(), out, 'fasta')
  tree.write(outfile=treeout, format=5)

if __name__ == '__main__':
  main("$reftree", "$refalignment", "${reftree.simpleName}.trim.tree", "${refalignment.simpleName}.trim.aln", "${refalignment.simpleName}.tax.map", "$remove_taxids", "$precompiled_taxa")
