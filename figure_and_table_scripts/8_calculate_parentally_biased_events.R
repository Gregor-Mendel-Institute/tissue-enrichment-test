prefix='u_'
sourcedir="L:/members/Schon/TPC_revisions/datasets_generated/"
genecounts=gsub(paste('^',prefix,'(GSM|SRX)[^_]+_',sep=''),prefix,gsub('\\.counts','',list.files(paste(sourcedir,prefix,"genecounts/",sep=""))))
guide_table <- read.delim(paste(sourcedir,"guide_table.txt",sep=''), stringsAsFactors=FALSE)
rownames(guide_table)=guide_table[,1]
rownames(guide_table)=paste(prefix,rownames(guide_table),sep='')
genecounts=genecounts[genecounts%in%rownames(guide_table)]


########################################################
# Fisher's Exact Test on all pairs of reciprocal cross #
########################################################
guide_table=guide_table[which(!is.na(guide_table$reciprocal.partner)),]
genecounts=genecounts[genecounts%in%rownames(guide_table)]

minimum_required_reads=5
minimum_fold=1
minimum_pval=.01

pvals_recip=list()
mat_recip=list()
pat_recip=list()

for(i in genecounts){
  ri=paste(prefix,guide_table[i,5],sep='')
  cross_sub=cbind(cross_table[[i]],cross_table[[ri]])
  cross_sub=cross_sub[which(rowSums(cross_sub)>=minimum_required_reads&rowSums(is.na(cross_sub))==0),]
  colnames(cross_sub)=paste(colnames(cross_sub),rep(1:2,each=2),sep="")

  mat1=grep(paste(guide_table[i,'maternal'],1,sep=""),colnames(cross_sub))
  mat2=grep(paste(guide_table[ri,'maternal'],2,sep=""),colnames(cross_sub))
  pat1=grep(paste(guide_table[i,'paternal'],1,sep=""),colnames(cross_sub))
  pat2=grep(paste(guide_table[ri,'paternal'],2,sep=""),colnames(cross_sub))

  if(guide_table[i,'expected.ratio']==2){
    cross_sub[,mat1]=round(cross_sub[,mat1]/2)
    cross_sub[,mat2]=round(cross_sub[,mat2]/2)
  }

  matfold=names(which(log2((cross_sub[,mat1]+1)/(cross_sub[,pat1]+1))>=minimum_fold & log2((cross_sub[,mat2]+1)/(cross_sub[,pat2]+1))>=minimum_fold))
  patfold=names(which(log2((cross_sub[,pat1]+1)/(cross_sub[,mat1]+1))>=minimum_fold & log2((cross_sub[,pat2]+1)/(cross_sub[,mat2]+1))>=minimum_fold))

  Xsq=p.adjust(apply(cross_sub,1,function(rw)fisher.test(rbind(rw[1:2],rw[3:4]))$p.value),method = 'fdr')
  pvals_recip[[i]]=Xsq

  signif=names(which(pvals_recip[[i]]<minimum_pval))
  signif_mat=signif[signif%in%matfold]
  signif_pat=signif[signif%in%patfold]

  mat_recip[[i]]=Xsq[signif_mat]
  pat_recip[[i]]=Xsq[signif_pat]
  cat(gsub("\\.counts","",i),length(signif),length(signif_mat),length(signif_pat),sep="\t")
  cat("\n")
  pdf(paste(sourcedir,'/fisher_plots/',gsub("\\.counts","",i),".pdf",sep=""),8,8)
  smoothScatter(x=log2((cross_sub[,mat1]+1)/(cross_sub[,pat1]+1)),y=log2((cross_sub[,mat2]+1)/(cross_sub[,pat2]+1)),nbin = 256,nrpoints = round(nrow(cross_sub)*.25),
                pch=20,col = 'gray',cex=.5,colramp = colorRampPalette(c('white','gray98','gray60','black')),las=1,
                xlim=c(-10,10),ylim=c(-10,10),
                ylab=paste("log2(",paste(as.character(guide_table[ri,c('maternal','paternal')]),collapse="/"),"),  ",gsub("\\.counts","",ri),sep=""),
                xlab=paste("log2(",paste(as.character(guide_table[i,c('maternal','paternal')]),collapse="/"),"),  ",gsub("\\.counts","",i),sep=""))
  title(main=gsub("\\.counts","",i))
  abline(h=c(-minimum_fold,minimum_fold),lty=2,col="gray")
  abline(v=c(-minimum_fold,minimum_fold),lty=2,col="gray")
  abline(h=0)
  abline(v=0)
  points(x=log2((cross_sub[signif_mat,mat1]+1)/(cross_sub[signif_mat,pat1]+1)),y=log2((cross_sub[signif_mat,mat2]+1)/(cross_sub[signif_mat,pat2]+1)),pch=20,col="indianred1",cex=1)
  points(x=log2((cross_sub[signif_pat,mat1]+1)/(cross_sub[signif_pat,pat1]+1)),y=log2((cross_sub[signif_pat,mat2]+1)/(cross_sub[signif_pat,pat2]+1)),pch=20,col="lightskyblue",cex=1)
  text(c(10,-10),c(10,-10),labels = c(length(signif_mat),length(signif_pat)),col=c("indianred1","lightskyblue"),cex = 1.5,font=2)
  dev.off()
}

embryo_crosses = c("u_2cell_ColxCvi",
                   "u_8cell_ColxCvi",
                   "u_eg_ColxCvi",
                   "u_8cell_pilot_ColxCvi",
                   "u_Col-WTxLer-WT_Embryo",
                   "u_gehring_ColLer_emb",
                   "u_ColxCvi_embryo",
                   "u_CvixLer_embryo")

endosperm_crosses = c("u_Col-WTxLer-WT_Endosperm",
                      "u_Col-WTxLer-WT_full-Endosperm",
                      "u_gehring_ColLer_endo",
                      "u_ColxCvi_endosperm_1",
                      "u_ColxCvi_endosperm_2",
                      "u_CvixLer_endosperm",
                      "u_ColxLer_endosperm_2",
                      "u_ColxLer_endosperm_4",
                      "u_ColxLer_endosperm_5",
                      "u_ColxCvi_endosperm_3",
                      "u_CvixLer_endosperm_2",
                      "u_CvixLer_endosperm_3")

mat_recip = mat_recip[grep('whole_seed',names(mat_recip),invert = T,value = T)]
pat_recip = pat_recip[grep('whole_seed',names(pat_recip),invert = T,value = T)]

allmaternal = sort(unique(unlist(lapply(mat_recip,names))))
allpaternal = sort(unique(unlist(lapply(pat_recip,names))))

allevents = sort(union(allpaternal,allmaternal))

emb_mat = table(unlist(lapply(mat_recip[embryo_crosses],names)))
emb_pat = table(unlist(lapply(pat_recip[embryo_crosses],names)))
endo_mat = table(unlist(lapply(mat_recip[endosperm_crosses],names)))
endo_pat = table(unlist(lapply(pat_recip[endosperm_crosses],names)))

sink('imprint_AGIS')
TAIRnames(allevents)
sink()
write.table(snps_per_gene[allevents,],'imprint_snps.tsv',quote=F,sep='\t')

for(i in embryo_crosses){
  event_hits=rep('',length(allevents))
  event_hits[allevents%in%names(mat_recip[[i]])]='m'
  event_hits[allevents%in%names(pat_recip[[i]])]='p'
  cat(i,'\n')
  writeClipboard(as.character(event_hits))
  readline()
}
for(i in endosperm_crosses){
  event_hits=rep('',length(allevents))
  event_hits[allevents%in%names(mat_recip[[i]])]='m'
  event_hits[allevents%in%names(pat_recip[[i]])]='p'
  cat(i,'\n')
  writeClipboard(as.character(event_hits))
  readline()
}

write.table(emb_pat,'emb_pat.tsv',quote=F,sep='\t')
write.table(emb_mat,'emb_mat.tsv',quote=F,sep='\t')
write.table(endo_pat,'endo_pat.tsv',quote=F,sep='\t')
write.table(endo_mat,'endo_mat.tsv',quote=F,sep='\t')

