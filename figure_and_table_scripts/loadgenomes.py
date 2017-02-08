from rc import rc
loadSA = False
loadLCP = False
directory = "L:/members/Schon/R/Arabidopsis Genome/"

FASTA = {}
SA = {}
LCP = {}

for i in ['1','2','3','4','5','C','M']:
    print("Loading FASTA for chr"+i+"...")
    currentfile=open(directory+"TAIR10_chr"+i+".fas").readlines()
    FASTA[i]="".join([j.rstrip() for j in currentfile[1:len(currentfile)]])
    if loadSA:
        print("Loading SA for chr"+i+"...")
        currentfile=open(directory+"Ath"+i+".SA").readlines()
        SA[i]=[int(j) for j in currentfile]
        assert len(SA[i]) == len(FASTA[i])
    if loadLCP:
        print("Loading LCP for chr"+i+"...")
        currentfile=open(directory+"Ath"+i+".LCP").readlines()
        SA[i]=[int(j) for j in currentfile]
        assert len(SA[i]) == len(FASTA[i])-1


def grabseq(chrom,starts,stops,strand):
    if type(starts) is int and type(stops) is int:
        x=eval("Ath"+str(chrom))[(starts-1):stops]
    else:
        assert len(starts)==len(stops)
        x=[]
        for strt,stp in tuple(zip(starts,stops)):
            x.append(eval("Ath"+str(chrom))[(strt-1):stp])
        x="".join(x)
    if strand=="-":x=rc(x)
    return x

