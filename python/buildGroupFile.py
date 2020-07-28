import argparse

parser = argparse.ArgumentParser(description='Generate a group file from VEP output.')
parser.add_argument('-i', '--input', required=True,
                    help='input variant file (VEP output, standard VEP format)')
parser.add_argument('-o', '--output', required=True,
                    help='output group file, will overwrite existing files')
args = parser.parse_args()

inFileName = args.input #'variant_effect_filtered.txt'
outFileName = args.output #'group_file.txt'

print("Collecting variants from: " + inFileName)
count = 0
variants_by_gene = {}
with open(inFileName) as fp:
    for line in fp:
        if line.startswith("#"):
            continue
        count += 1
        fields = line.split("\t")
        gene = fields[3]
        variant = fields[0]
	variant = variant.replace("_", ":").replace("/", ":")
        if gene in variants_by_gene:
            gene_variants = variants_by_gene[gene]
        else:
            gene_variants = []
            variants_by_gene[gene] = gene_variants
        gene_variants.append(variant)

print("Got " + str(count) + " variants in " +
      str(len(variants_by_gene)) + " genes.")

outFile = open(outFileName, "w")
for gene in variants_by_gene:
    variants = variants_by_gene[gene]
    outFile.write(gene + "\t" + "\t".join(variants) + "\n")
outFile.close()

print("Wrote group file: " + outFileName)



