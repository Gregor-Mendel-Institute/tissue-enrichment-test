def rc(mir):
    mir=mir.upper()
    rc=["0"]*len(mir)
    for i in range(len(mir)-1,-1,-1):
        if(mir[len(mir)-1-i]=="A"):rc[i]="T"
        elif(mir[len(mir)-1-i]=="C"):rc[i]="G"
        elif(mir[len(mir)-1-i]=="G"):rc[i]="C"
        elif(mir[len(mir)-1-i]=="T" or mir[len(mir)-1-i]=="U"):rc[i]="A"
        else:rc[i]="."
    return("".join(rc))
