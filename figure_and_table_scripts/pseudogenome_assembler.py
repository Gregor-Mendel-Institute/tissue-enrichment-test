import os,re
from loadgenomes import *
from table import *
from rc import *
ecotype_ID={}
ecotype_ID['6909']='Col'
ecotype_ID['6911']='Cvi'
ecotype_ID['7213']='Ler'

dataset_folder='../datasets_external/'
onlySNPs=False
GTFfile='Araport11_GFF3_genes_transposons.201606.gtf'

if onlySNPs==False:
    print("###\nWARNING: Currently only supports pseudogenome construction with SNPs.")
    print("Changing options to only use SNPs.\n###")
    onlySNPs=True

pseudo={}
polymorphisms={}
indel_index={}

for k,v in ecotype_ID.items():
    intersect=open(dataset_folder+'intersection_'+k+'_'+v+'.vcf')
    polymorphisms[v]={}
    curchrom=0
    print("Loading polymorphism index for "+v+"...")
    for line in intersect:
        if line[0] in ['#','C']:continue
        l=line.rstrip().split('\t')
        chrom,pos,ID,ref,alt=l[0:5]
        if chrom!=curchrom:
            polymorphisms[v][chrom]={}
            curchrom=chrom
        if onlySNPs:
            if len(ref)==1 and len(alt)==1:
                polymorphisms[v][chrom][pos]=(ref,alt)
        else:
            polymorphisms[v][chrom][pos]=(ref,alt)
    lens=[len(polymorphisms[v][c]) for c in ['1','2','3','4','5']]
    print("Number of polymorphisms:")
    for i in ['1','2','3','4','5']:
        print('chr'+i+'\t'+str(lens[int(i)-1]))
        
            
print("Workspace Loaded.\n#######################################\nMaking pseudogenomes...")
for k,v in ecotype_ID.items():
    pseudo[v]={}
    print("Assembing pseudogenome for "+v+"...")
    for chrom in ['1','2','3','4','5']:
        print('\tchr'+chrom)
        pseudo[v][chrom]=list(FASTA[chrom])
        for pk,pv in polymorphisms[v][chrom].items():
            old,new=pv
            assert pseudo[v][chrom][int(pk)-1]==old
            pseudo[v][chrom][int(pk)-1]=new
    print("Writing FASTA file for "+v+"...")
    outfile=open("pseudo_"+v+".fa","w",newline='')
    #for chrom in ['1','2','3','4','5']:
    #    outfile.write(">"+v+"_chr"+chrom+"\n")
    #    outfile.write("".join(pseudo[v][chrom])+"\n")
    #for chrom in ['C','M']:
    #    outfile.write(">"+v+"_chr"+chrom.lower()+"\n")
    #    outfile.write(FASTA[chrom]+"\n")
    outfile.close()
    print("Pseudogenome saved to pseudo_"+v+".fa")
    print("Making "+v+"-specific version of "+GTFfile)
    GTF=open(dataset_folder+GTFfile)
    outGTF=open("Araport11_"+v+".gtf","w",newline='')
    #for line in GTF:
    #    l=line.rstrip().split('\t')
    #    l[0]="_".join([v,l[0].lower()])
    #    if not l[8].endswith(';'):
    #        l[8]=l[8]+';'
    #    m=re.search('(transcript_id \".+)(\..+?\"; gene_id \".+?)(\";)',l[8])
    #    l[8]=m.groups()[0]+v.lower()+m.groups()[1]+v.lower()+m.groups()[2]
    #    outGTF.write("\t".join(l)+'\n')
    outGTF.close()
    print("Ecotype GTF saved to Araport11_"+v+".gtf")
print("#######################################")
for v1,v2 in [('Col','Ler'),('Col','Cvi'),('Cvi','Ler')]:
    print("Making hybrid pseudogenome for "+v1+' and '+v2+'...')
    for fileheader,filetype in [('Araport11_','.gtf'),('pseudo_','.fa')]:
        file1=open(fileheader+v1+filetype).read()
        file2=open(fileheader+v2+filetype).read()
        hybrid=open(fileheader+'_'.join([v1,v2])+filetype,'w',newline='')
        hybrid.write(file1+file2)
        hybrid.close()
        print('\tSaved '+fileheader+'_'.join([v1,v2])+filetype)

for combos in [('6909','7213'),('6909','6911'),('7213','6911')]:
    print('Total SNPs between'+ecotype_ID[combos[0]]+' and '+ecotype_ID[combos[1]]+':')
    different_locations = 0
    for chrom in ['1','2','3','4','5']:
        different_locations += sum([i[0]!=i[1] for i in zip(pseudo[ecotype_ID[combos[0]]][chrom],pseudo[ecotype_ID[combos[1]]][chrom])])
    print(str(different_locations))
