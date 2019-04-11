import glob

all_tax = {}
for tax_map in glob.glob("references_trimmed/*.tax.map"):
     for line in open(tax_map):
         line = line.strip().split()
         all_tax[line[0].split('.')[0]] = line[1]

with open('precompiled.tax.map', 'w') as out:
    for tax, lin in all_tax.items():
        print("{}\t{}".format(tax, lin), file=out)
