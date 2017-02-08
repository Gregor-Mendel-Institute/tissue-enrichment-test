from loadgenomes_single_fasta import *
from rc import *
print(FASTA.keys())
import os,sys,re
def which(x,value=True):
    return [a for a,b in enumerate(x) if b==value]
def notwhich(x,value=0):
    return [a for a,b in enumerate(x) if b!=value]
indexing = 1
samfile=open(sys.argv[1])
#SNPfile=open(sys.argv[2])
#samfile=open('../LerCol_24cell_11.sam')
ler_out=open(sys.argv[1]+'_ler','w')
col_out=open(sys.argv[1]+'_col','w')
SNPfile=open('../autran_replication/Ler_SNPs_Ecker.txt')
SNPlocs={}
for chrom in ['1','2','3','4','5','c','m']:
    SNPlocs[chrom]={}
SNPfile.readline()
matches_tair8=[0,0]
for line in SNPfile:
    eco,chrom,loc,old,new,snptype=line.split('\t')[0:6]
    if FASTA[chrom[-1]][int(loc)-1]==old:
        SNPlocs[chrom[-1]][str(int(loc)-1)]=(old,new)
        matches_tair8[0]+=1
    else:
        matches_tair8[1]+=1
print('match: ',matches_tair8[0],' nonmatch: ',matches_tair8[1])
MM_col=16
ignoreIDN=True
SNPtypes={}
assertion_counter=0
for line in samfile:
    if line[0]=='@':continue
    l=line.rstrip().split('\t')
    name,flag,chrom,pos,score,cigar,a,b,c,seq,qual,NH,HI,AS,nM=l
    if chrom[-1] not in SNPlocs.keys():
        continue
    if flag=='0':
        strand='plus'
    elif flag=='16':
        strand='minus'
    else:
        continue
    start_pos=int(pos)-1
    end_pos=start_pos+len(seq)
    cigar_nums=re.split('[A-Z]+',cigar)
    cigar_nums=[int(i) for i in cigar_nums[0:len(cigar_nums)-1]]
    cigar_chars=re.split('[0-9]+',cigar)
    cigar_chars=cigar_chars[1:len(cigar_chars)]
    nucleotides=list(range(start_pos,end_pos))
    if ignoreIDN:
        if 'I' in cigar_chars or 'D' in cigar_chars or 'N' in cigar_chars:
            continue
    if 'I' in cigar_chars:
        end_pos += -sum([cigar_nums[a] for a in which(cigar_chars,'I')])
        nucleotides###        
    if 'N' in cigar_chars:
        end_pos += sum([cigar_nums[a] for a in which(cigar_chars,'N')])
        nucleotides###
    if 'D' in cigar_chars:
        end_pos += sum([cigar_nums[a] for a in which(cigar_chars,'I')])
        
    is_a_snp=[str(i) in SNPlocs[chrom[-1]] for i in nucleotides]
    #print(' '.join([chrom[-1],str(start_pos),str(end_pos)]))
    if any(is_a_snp):
        
        #print(line)
        #print([' '.join([FASTA[chrom[-1]][start_pos+i-2:start_pos+i+3],str(SNPlocs[chrom[-1]][str(start_pos+i)]),seq[i]]) for i in which(is_a_snp)],sep='\n')
        #print(SNPlocs[chrom[-1]][str(start_pos+which(is_a_snp)[0])])
        assert assertion_counter<=100
        if all([SNPlocs[chrom[-1]][str(start_pos+i)][0]==seq[i] for i in which(is_a_snp)]):
            readtype='Col'
            SNPtypes['Col']=SNPtypes.get('Col',0)+1
            col_out.write(line)
        elif all([SNPlocs[chrom[-1]][str(start_pos+i)][1]==seq[i] for i in which(is_a_snp)]):
            readtype='Ler'
            SNPtypes['Ler']=SNPtypes.get('Ler',0)+1
            ler_out.write(line)
        else:
            readtype='unresolved'
            SNPtypes['unresolved']=SNPtypes.get('unresolved',0)+1
        #print(seq,readtype,' '.join([seq[i] for i in which(is_a_snp)]))
print(SNPtypes)
col_out.close()
ler_out.close()
