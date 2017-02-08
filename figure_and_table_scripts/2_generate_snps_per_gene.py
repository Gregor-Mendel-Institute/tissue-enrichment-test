import re,os
from FASTA_tools import *

Col = import_FASTA('../datasets_generated/pseudo_Col.fa')
Ler = import_FASTA('../datasets_generated/pseudo_Ler.fa')
Cvi = import_FASTA('../datasets_generated/pseudo_Cvi.fa')

annotation_file = '../datasets_external/Araport11_GFF3_genes_transposons.201606.gtf'
genes = {}

output_file = open('../datasets_generated/snps_per_gene.tsv','w')

for line in open(annotation_file):
    l = line.rstrip().split('\t')
    if l[2] != 'exon':
        continue
    chr_number = l[0][-1]
    if chr_number not in ['1','2','3','4','5']:
        continue
    gene = re.search('gene_id \"([^\"]+)\"',l[-1]).groups()[0]
    genes[gene] = genes.get(gene,[])+list(range(int(l[3]),int(l[4])+1))

gene_names = sorted(list(genes.keys()))

output_file.write('gene\tCol_Ler\tCol_Cvi\tCvi_Ler\n')
for g in gene_names:
    genes[g] = list(sorted(set(genes[g])))
    chr_number = g[2]
    Col_seq = ''.join([Col['Col_chr'+chr_number][j] for j in genes[g]])
    Ler_seq = ''.join([Ler['Ler_chr'+chr_number][j] for j in genes[g]])
    Cvi_seq = ''.join([Cvi['Cvi_chr'+chr_number][j] for j in genes[g]])
    Col_Ler = sum([i[0]!=i[1] for i in zip(Col_seq,Ler_seq)])
    Col_Cvi = sum([i[0]!=i[1] for i in zip(Col_seq,Cvi_seq)])
    Cvi_Ler = sum([i[0]!=i[1] for i in zip(Cvi_seq,Ler_seq)])
    output_file.write('\t'.join([g,str(Col_Ler),str(Col_Cvi),str(Cvi_Ler)])+'\n')
output_file.close()    
