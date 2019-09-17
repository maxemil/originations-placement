#!/usr/bin/env nextflow

params.fasta = ""
params.reference_alignments = ""
params.reference_trees = ""
params.remove_taxids = ""
params.og_hierarchy = ""
params.precompiled_taxa = ""

params.placment_dir = "placements"
params.references_trimmed_dir = "references_trimmed"
params.queries_dir = "halo_queries"

remove_taxids = Channel.from(file(params.remove_taxids))
og_hierarchy = Channel.from(file(params.og_hierarchy))
precompiled_taxa = Channel.from(file(params.precompiled_taxa))
fasta = Channel.fromPath(params.fasta)
reference_trees = Channel.fromPath(params.reference_trees).map{it -> tuple(it.simpleName, it)}

Channel.fromPath(params.reference_alignments).into{reference_alignments_queries; reference_alignments_tree}

references = reference_alignments_queries.map{it -> tuple(it.simpleName, it)}.join(reference_trees)

// process filterFastaRecs {
//   input:
//   file fasta from fasta
//   file taxa_selection from select_taxids.first()
//
//   output:
//   file "${fasta.simpleName}.clean.fasta" into cleaned_fasta
//
//   script:
//   """
//   #! /opt/conda/envs/placement/bin/python3
//   from Bio import SeqIO
//
//   taxa_selection = [t.strip() for t in open("$taxa_selection")]
//   with open("${fasta.simpleName}.clean.fasta", 'w') as out:
//       for rec in SeqIO.parse("$fasta", 'fasta'):
//           if not any([rec.id.startswith(t) for t in taxa_selection]):
//               print("remove seq {}".format(rec.id))
//           else:
//               SeqIO.write(rec, out, 'fasta')
//   """
// }

process matchRefQuery {
  input:
  file og_hierarchy from og_hierarchy.first()
  file fasta from fasta

  output:
  file "${fasta.simpleName}.*.fasta" into matched_fasta mode flatten
  file "*.og" into reference_ogs mode flatten
  tag {"${fasta.simpleName}"}

  script:
  """
  for og in \$(grep ${fasta.simpleName} $og_hierarchy | cut -f 1);
  do
    cp $fasta ${fasta.simpleName}.\$og.fasta
    touch \$og.og
  done
  """
}

unique_og = reference_ogs.unique{it.simpleName}

process resolvePolytomies {
  input:
  set val(id), file(refalignment), file(reftree) from references
  each og from unique_og
  file remove_taxids from remove_taxids.first()
  file precompiled_taxa from precompiled_taxa.first()

  output:
  set file("${reftree.simpleName}.trim.tree"), file("${refalignment.simpleName}.trim.aln") into binary_reftrees optional true
  file "${refalignment.simpleName}.trim.aln" into trim_refaln optional true
  file "${refalignment.simpleName}.tax.map" into tax_map optional true
  file "references_summary.txt" into ref_summary

  tag {"${reftree.simpleName}"}
  publishDir "$params.references_trimmed_dir", mode: 'copy'

  when:
  "$id" == "${og.simpleName}"

  script:
  template 'prepare_eggnog_ref_OG.py'
}
ref_summary.collectFile(name: 'references_summary.txt', storeDir: "$workflow.launchDir")

process addQueries {
  input:
  each fasta from matched_fasta
  file refalignment from trim_refaln

  output:
  file("${fasta.baseName}.aln") into alignments_placement

  tag {"${fasta.baseName}"}
  cpus 4
  publishDir "$params.queries_dir", mode: 'copy'

  when:
  "${refalignment.simpleName}" == "${fasta.baseName.tokenize('.')[1]}"

  script:
  """
  trimal -in $refalignment -out ${refalignment.baseName}.90.aln -gt 0.1 -fasta -keepseqs
  mafft --addfragments $fasta --thread ${task.cpus} --keeplength ${refalignment.baseName}.90.aln > ${fasta.baseName}.refquery.aln
  trimal -in ${refalignment.baseName}.90.aln -out ${refalignment.simpleName}.phy -phylip
  trimal -in ${fasta.baseName}.refquery.aln -out ${fasta.baseName}.refquery.phy -phylip
  epa-ng --split ${refalignment.simpleName}.phy ${fasta.baseName}.refquery.phy
  mv query.fasta ${fasta.baseName}.aln
  """
}

process evaluateTree {
  input:
  set file(tree), file(refmsa) from binary_reftrees

  output:
  set file("$tree"), file("${refmsa.baseName}.90.aln.raxml.bestModel"), file("${refmsa.baseName}.90.aln") into tree_model

  tag {"${tree.simpleName}"}
  cpus 4
  publishDir "$params.references_trimmed_dir", pattern: '*.bestModel', mode: 'copy'

  script:
  """
  trimal -in $refmsa -out ${refmsa.baseName}.90.aln -gt 0.1 -fasta -keepseqs
  raxml-ng --evaluate --tree $tree --opt-branches off --msa ${refmsa.baseName}.90.aln --model LG+F+I --threads ${task.cpus} --force
  """
}

process placeQueries {
  input:
  set file(tree), file(modelinfo), file(refmsa) from tree_model
  each querymsa from alignments_placement

  output:
  file "${querymsa.baseName}.jplace" into jplace_files

  tag {"${querymsa.simpleName} - ${tree.simpleName}"}
  cpus 4
  publishDir "$params.placment_dir", mode: 'copy'

  when:
  "${tree.simpleName}" == "${querymsa.baseName.tokenize('.')[1]}"

  script:
  """
  epa-ng --ref-msa $refmsa --tree $tree --query $querymsa --threads ${task.cpus} --verbose --model $modelinfo
  mv epa_result.jplace ${querymsa.baseName}.jplace
  """
}


process assignTaxonomy {
  input:
  file jplace from jplace_files
  each tax_map from tax_map

  output:
  file "${jplace.baseName}.csv" into profiles
  file "${jplace.baseName}.assign" into assignments
  file "${jplace.baseName}.svg" into svg_trees
  file "${jplace.baseName}.newick" into newick_trees

  tag {"${jplace.simpleName} - ${tax_map.simpleName}"}
  publishDir "$params.placment_dir", mode: 'copy'

  when:
  "${tax_map.simpleName}" == "${jplace.baseName.tokenize('.')[1]}"

  script:
  """
  gappa analyze assign --jplace-path $jplace --taxon-file $tax_map --threads ${task.cpus}
  mv profile.tsv ${jplace.baseName}.csv
  mv per_query.tsv ${jplace.baseName}.assign
  gappa analyze visualize-color --jplace-path $jplace --write-svg-tree --threads ${task.cpus}
  mv tree.svg ${jplace.baseName}.svg
  gappa analyze graft --name-prefix 'Q_' --jplace-path $jplace --threads ${task.cpus}
  """
}
