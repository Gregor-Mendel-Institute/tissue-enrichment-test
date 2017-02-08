def table(x):
    y=sorted(list(set(x)))
    z={}
    for i in y:
        z[i]=len([1 for b in x if b==i])
    return sorted([(v,k) for k,v in z.items()],reverse=True)
