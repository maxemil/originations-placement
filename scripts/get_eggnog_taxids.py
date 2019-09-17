from ete3 import ncbi_taxonomy, Tree
ncbi = ncbi_taxonomy.NCBITaxa()

methanotecta = ['Archaeoglobi', 'Methanonatronarchaeia', 'Halobacteria', 'Methanomicrobia', 'Marine Group IV']
with open('Methanotecta.taxid', 'w') as out:
    for line in open('eggnog4.species_list.txt'):
        if not line.startswith('#'):
            lineage = ncbi.get_lineage(int(line.split('\t')[1]))
            for l in lineage:
                if ncbi.get_taxid_translator([l])[l] in methanotecta:
                    print(line.split('\t')[1], file=out)
