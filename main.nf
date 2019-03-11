#!/usr/bin/env nextflow

params.fasta = ""
params.reference_alignments = ""
params.reference_trees = ""
params.taxid = ""
params.og_hierarchy = ""

taxids = Channel.from(file(params.taxid))
og_hierarchy = Channel.from(file(params.og_hierarchy))
fasta = Channel.fromPath(params.fasta)
reference_trees = Channel.fromPath(params.reference_trees)
Channel.fromPath(params.reference_alignments).into{reference_alignments_queries; reference_alignments_tree}


process matchRefQuery {
  input:
  file og_hierarchy from og_hierarchy.first()
  file fasta from fasta

  output:
  file "${fasta.simpleName}.*.fasta" into matched_fasta mode flatten

  tag {"${fasta.simpleName}"}

  script:
  """
  for og in \$(grep ${fasta.simpleName} $og_hierarchy | cut -d'@' -f 1);
  do
    cp $fasta ${fasta.simpleName}.\$og.fasta
  done
  """
}


process resolvePolytomies {
  input:
  file reftree from reference_trees
  each refalignment from reference_alignments_queries
  file remove_taxids from taxids.first()

  output:
  set file("${reftree.simpleName}.trim.tree"), file("${refalignment.simpleName}.trim.aln") into binary_reftrees
  file "${refalignment.simpleName}.trim.aln" into trim_refaln
  file "${refalignment.simpleName}.tax.map" into tax_map

  tag {"${reftree.simpleName}"}

  when:
  "${reftree.simpleName}" == "${refalignment.simpleName}"

  script:
  template 'prepare_eggnog_ref_OG.py'
}

process addQueries {
  input:
  each fasta from matched_fasta
  file refalignment from trim_refaln

  output:
  file("${fasta.baseName}.aln") into alignments_placement

  tag {"${fasta.baseName}"}

  when:
  "${refalignment.simpleName}" == "${fasta.baseName.tokenize('.')[1]}"

  script:
  """
  mafft --add $fasta --thread 20 --keeplength $refalignment > ${fasta.baseName}.refquery.aln
  trimal -in $refalignment -out ${refalignment.simpleName}.phy -phylip
  trimal -in ${fasta.baseName}.refquery.aln -out ${fasta.baseName}.refquery.phy -phylip
  epa-ng --split ${refalignment.simpleName}.phy ${fasta.baseName}.refquery.phy
  mv query.fasta ${fasta.baseName}.aln
  """
}

process evaluateTree {
  input:
  set file(tree), file(refmsa) from binary_reftrees

  output:
  set file("$tree"), file("${refmsa}.raxml.bestModel"), file("$refmsa") into tree_model

  tag {"${tree.simpleName}"}

  script:
  """
  raxml-ng --evaluate --tree $tree --opt-branches off --msa $refmsa --model LG+F+I --threads 3
  """
}

process placeQueries {
  input:
  set file(tree), file(modelinfo), file(refmsa) from tree_model
  each querymsa from alignments_placement

  output:
  file "${querymsa.baseName}.jplace" into jplace_files

  tag {"${querymsa.simpleName} - ${tree.simpleName}"}

  when:
  "${tree.simpleName}" == "${querymsa.baseName.tokenize('.')[1]}"

  script:
  """
  epa-ng --ref-msa $refmsa --tree $tree --query $querymsa --threads 40 --verbose --model $modelinfo
  mv epa_result.jplace ${querymsa.baseName}.jplace
  """
}


process assignTaxonomy {
  input:
  file jplace from jplace_files
  each tax_map from tax_map

  output:

  tag {"${jplace.simpleName} - ${tax_map.simpleName}"}

  when:
  "${tax_map.simpleName}" == "${jplace.baseName.tokenize('.')[1]}"

  script:
  """
  gappa analyze assign --jplace-path $jplace --taxon-file $tax_map --threads ${task.cpus}
  mv profile.csv ${jplace.simpleName}.csv
  mv per_pquery_assign ${jplace.simpleName}.assign
  gappa analyze visualize-color --jplace-path $jplace --write-svg-tree --threads ${task.cpus}
  mv tree.svg ${jplace.simpleName}.svg
  gappa analyze graft --name-prefix 'Q_' --jplace-path $jplace --threads ${task.cpus}
  """
}
