import re
FASTA_file='L:/members/Schon/lab_files/genomes/TAIR10_ERCC.fas'
out_file="L:/members/Schon/genome_matches.bed"

def import_FASTA(FASTA_file):
    genome={}
    chromosomes={}
    chrom='none'
    genome_file=open(FASTA_file)
    for line in genome_file:
        line=line.rstrip()
        if line[0]=='>':
            if chrom!='none':
                genome[chrom]=''.join(x)
            chrom=line[1:len(line)]
            x=[]
            continue
        x.append(line)
    genome[chrom]=''.join(x)
    for k,v in genome.items():
        # assert len(v) == chromosomes[k]
        chromosomes[k]=len(v)
    return genome

IUPAChash={}
IUPACcomp={}
values=['A','C','G','T','T',
      '[AC]','[AG]','[AT]','[CG]','[CT]','[GT]',
      '[ACG]','[ACT]','[AGT]','[CGT]','.']
keys=['A','C','G','T','U',
        'M','R','W','S','Y','K',
        'V','H','D','B','N']
complements=['T','G','C','A','A',
             'K','R','W','S','Y','M',
             'B','D','H','V','N']
for k,v,c in zip(keys,values,complements):
     IUPAChash[k]=v
     IUPACcomp[k]=c

def rc(sequence):
    return ''.join(reversed([IUPACcomp[i] for i in sequence.upper()]))

def to_regex(sequence):
    return ''.join([IUPAChash[i] for i in sequence.upper()])

def genome_matches(sequence):
    out_bed=open(out_file,'w')
    seq1=to_regex(sequence)
    seq2=to_regex(rc(sequence))
    print(seq1,seq2)
    hit_number=1
    for chrom in ['Ath_chr1','Ath_chr2','Ath_chr3','Ath_chr4','Ath_chr5','Ath_chrc','Ath_chrm']:
        hits_sense=[i.span() for i in re.finditer(seq1,genome[chrom])]
        for h in hits_sense:
            out_bed.write('\t'.join([chrom,str(h[0]),str(h[1]),'match'+str(hit_number),'0','+'])+'\n')
            hit_number+=1
        hits_antisense=[i.span() for i in re.finditer(seq2,genome[chrom])]
        for h in hits_sense:
            out_bed.write('\t'.join([chrom,str(h[0]),str(h[1]),'match'+str(hit_number),'0','-'])+'\n')
            hit_number+=1
    out_bed.close()
