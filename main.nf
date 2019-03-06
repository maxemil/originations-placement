#!/usr/bin/env nextflow

params.fasta = ""
params.reference_alignments = ""
params.reference_trees = ""
params.taxid = ""

taxids = Channel.from(file(params.taxid))
fasta = Channel.fromPath(params.fasta)
reference_trees = Channel.fromPath(params.reference_trees)
Channel.fromPath(params.reference_alignments).into{reference_alignments_queries; reference_alignments_tree}

process resolvePolytomies {
  input:
  file reftree from reference_trees
  file refalignment from reference_alignments_queries
  file remove_taxids from taxids
  
  output:
  set file("${reftree.simpleName}.trim.tree"), file("${refalignment.simpleName}.trim.aln") into binary_reftrees
  file "${refalignment.simpleName}.trim.aln" into trim_refaln
  file "${refalignment.simpleName}.tax.map" into tax_map

  tag {"${reftree.simpleName}"}

  script:
  template 'prepare_eggnog_ref_OG.py'
}

process addQueries {
  input:
  file fasta from fasta
  file refalignment from trim_refaln

  output:
  file("${fasta.simpleName}.aln") into alignments_placement

  tag {"${fasta.simpleName}"}

  script:
  """
  mafft --add $fasta --thread 20 --keeplength $refalignment > ${fasta.simpleName}.${refalignment.simpleName}.aln
  trimal -in $refalignment -out ${refalignment.simpleName}.phy -phylip
  trimal -in ${fasta.simpleName}.${refalignment.simpleName}.aln -out ${fasta.simpleName}.${refalignment.simpleName}.phy -phylip
  singularity exec -B /local:/local /local/one/people/MaxEmil/MarineGroupIV/23_placing_acquisitions/epa-ng epa-ng --split ${refalignment.simpleName}.phy ${fasta.simpleName}.${refalignment.simpleName}.phy
  mv query.fasta ${fasta.simpleName}.aln
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
  singularity exec -B /local:/local /local/one/people/MaxEmil/MarineGroupIV/23_placing_acquisitions/raxml-ng.simg raxml-ng --evaluate --tree $tree --opt-branches off --msa $refmsa --model LG+F+I --threads 3
  """
}

process placeQueries {
  input:
  set file(tree), file(modelinfo), file(refmsa) from tree_model
  file(querymsa) from alignments_placement

  output:
  file "${querymsa.simpleName}.jplace" into jplace_files

  tag {"${querymsa.simpleName}"}

  script:
  """
  singularity exec -B /local:/local /local/one/people/MaxEmil/MarineGroupIV/23_placing_acquisitions/epa-ng epa-ng --ref-msa $refmsa --tree $tree --query $querymsa --threads 40 --verbose --model $modelinfo
  mv epa_result.jplace ${querymsa.simpleName}.jplace
  """
}


process assignTaxonomy {
  input:
  file jplace from jplace_files
  file tax_map from tax_map

  output:

  tag {"${jplace.simpleName}"}

  script:
  """
  singularity exec -B /local:/local /local/one/people/MaxEmil/MarineGroupIV/23_placing_acquisitions/gappa gappa analyze assign --jplace-path $jplace --taxon-file $tax_map --threads ${task.cpus}
  mv profile.csv ${jplace.simpleName}.csv
  mv per_pquery_assign ${jplace.simpleName}.assign
  singularity exec -B /local:/local /local/one/people/MaxEmil/MarineGroupIV/23_placing_acquisitions/gappa gappa analyze visualize-color --jplace-path $jplace --write-svg-tree --threads ${task.cpus}
  mv tree.svg ${jplace.simpleName}.svg
  singularity exec -B /local:/local /local/one/people/MaxEmil/MarineGroupIV/23_placing_acquisitions/gappa gappa analyze graft --name-prefix 'Q_' --jplace-path $jplace --threads ${task.cpus}
  """
}
