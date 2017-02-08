import sys
filename=sys.argv[1]
outname=sys.argv[2]
f=open(filename)

reads={}
outfile=open(outname,'wb')
firstline=True
curlinetype=''
ISFASTA=False
for line in f:
    if line[0]=='@':
        seqID=line
        curlinetype='IDseq'
    elif line[0]=='+':
        IDqual=line
        qualID=line
        curlinetype='IDqual'
    elif curlinetype=='IDseq':
        seq=line
        reads[seq]=reads.get(seq,0)+1
        if reads[seq]==1:
            keepit=True
        else:
            keepit=False
    elif curlinetype=='IDqual':
        if len(line)==1:continue
        qual=line
        if keepit:
            outfile.write(seqID+seq+qualID+qual)
    elif firstline and line[0]==">":
        print 'Actually a FASTA file'
        ISFASTA=True
        break
    firstline=False

if ISFASTA:
    f.seek(0)
    for line in f:
        if line[0]=='>':
            seqID=line
            curlinetype='IDseq'
        elif curlinetype=='IDseq':
            seq=line
            reads[seq]=reads.get(seq,0)+1
            if reads[seq]==1:
                outfile.write(seqID+seq)

outfile.close()

#readslist=sorted([(v,k) for k,v in reads.items() if v>1],reverse=True)
#logfile=open(outname.replace('.fastq','')+'.mostabund','wb')
#for v,k in readslist:
#    logfile.write(k.rstrip()+'\t'+str(v)+'\n')
#logfile.close()

