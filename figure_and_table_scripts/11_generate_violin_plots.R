setwd('C:/Users/schon.admin/Desktop/TPC_revisions/datasets_generated/')

TPM = read.table('rsem_TPM.tsv',sep='\t',header = T,stringsAsFactors = F,row.names = 1)
colnames(TPM)=gsub('^X','',gsub('GaD_dcl1','dcl1_globular',gsub('GaD_wt','wt_globular',gsub('.','-',gsub('^(GSM|SRX)[^_]+_','',colnames(TPM)),fixed = T))))

belmonte_LCM=read.table('reference_atlas.tsv',header = T,row.names = 1)
reference_names=read.table('reference_atlas.description')
data_table=read.table('all_samples_in_study.tsv',header = T,stringsAsFactors = F,comment.char = '#',sep='\t',quote = '')
for(i in 1:ncol(belmonte_LCM)){
  colnames(belmonte_LCM)[i]=data_table[which(data_table[,"Alternate.Name"]==colnames(belmonte_LCM)[i]),1]
}
TPM = read.table('rsem_TPM.tsv',sep='\t',header = T,stringsAsFactors = F,row.names = 1)
colnames(TPM)=gsub('^X','',gsub('GaD_dcl1','dcl1_globular',gsub('GaD_wt','wt_globular',gsub('.','-',gsub('^(GSM|SRX)[^_]+_','',colnames(TPM)),fixed = T))))
belmonte_normalized=matrix(nrow=nrow(TPM),ncol=ncol(belmonte_LCM))
rownames(belmonte_normalized)=rownames(TPM)
colnames(belmonte_normalized)=colnames(belmonte_LCM)
in_new=rownames(belmonte_LCM)%in%rownames(belmonte_normalized)
for(i in colnames(belmonte_LCM)){
  belmonte_normalized[rownames(belmonte_LCM)[in_new],i]=belmonte_LCM[in_new,i]
}

allvalues=cbind(TPM,belmonte_normalized)
rm(belmonte_normalized)
for(i in 1:ncol(allvalues)){
  allvalues[,i]=allvalues[,i]/sum(allvalues[,i],na.rm = T)
}

x=read.table('emb_pat.tsv',header = T,stringsAsFactors = F)
emb_pat=x[,2]
names(emb_pat)=x[,1]
x=read.table('emb_mat.tsv',header = T,stringsAsFactors = F)
emb_mat=x[,2]
names(emb_mat)=x[,1]
x=read.table('endo_pat.tsv',header = T,stringsAsFactors = F)
endo_pat=x[,2]
names(endo_pat)=x[,1]
x=read.table('endo_mat.tsv',header = T,stringsAsFactors = F)
endo_mat=x[,2]
names(endo_mat)=x[,1]


contam_ratio=log2((rowSums(data_table[,18:19],na.rm = T)+1)/(rowSums(data_table[,c(15:17)],na.rm = T)+1))
names(contam_ratio)=data_table[,1]
smpls=which(data_table$Tissue.Type %in% c('endosperm','whole seed','seed coat','general seed coat','chalazal seed coat','peripheral endosperm','chalazal endosperm','micropylar endosperm') &
  data_table$Corresponding.Developmental.Stage=='linear cotyledon')
pdf('contam_ratio_barplot.pdf',useDingbats = F)
par(mar=c(8,4,2,2))
barplot(unlist(sort(contam_ratio[smpls])),las=2,cex.names=.6,space=0)
dev.off()
contamination_values=matrix(nrow=nrow(allvalues),ncol=length(smpls))
rownames(contamination_values)=rownames(allvalues)
colnames(contamination_values)=data_table[smpls,1]
for(i in 1:length(smpls)){
  contamination_values[,data_table[smpls[i],1]]=allvalues[,data_table[smpls[i],1]]
}
contam_ratio=contam_ratio[names(contam_ratio[smpls])[names(contam_ratio[smpls])%in%colnames(contamination_values)]]
contamination_values=contamination_values[,names(sort(contam_ratio))]
cr=as.numeric(sort(contam_ratio))
contam_cors=apply(contamination_values,1,function(x){
  if(length(which(is.na(x)))>=(length(x)-3)){NA}else{
    use_only=which(!is.na(x))
    cor(y=(x[use_only]),x=cr[use_only],method = 'spearman')}
})
contam_pvals=apply(contamination_values,1,function(x){
  if(length(which(is.na(x)))>=(length(x)-3)){NA}else{
    use_only=which(!is.na(x))
    cor.test(y=(x[use_only]),x=cr[use_only],method = 'spearman')$p.value}
})

########################################################################################
# OUR ANALYSIS - POTENTIAL IMPRINTED GENES
########################################################################################
allendomat=endo_mat
allendopat=endo_pat
dot=20
library(vioplot)
contam_mat=list()
contam_mat[['0']]=names(contam_cors[!names(contam_cors)%in%names(allendomat)])
for(i in c('1','2:5','6:12')){
  contam_mat[[as.character(i)]]=names(allendomat[which(allendomat%in%eval(parse(text=i)))])
}

pdf("maternal_imprinting.pdf",7,7,useDingbats = F)
realmaternals=list()
realmaternals[['2:5']]=c(agi("FWA"),agi("FIS2"))
realmaternals[['6:12']]=c(agi("^SDC$"),agi("MRU1"),agi("HDG9"))
scmaternals=list()
scmaternals[['2:5']]=c(agi("SHP1"),agi("SHP2"),agi("ELA1"),agi("^TT1$"),agi("^TT5$"),agi("^TT16$"))
scmaternals[['6:12']]=c(agi("STK"),agi("BAN"),agi("TT10"),agi("^TT3$"),agi("^TT19$"),agi("^TT4$"),agi("^TT7$"),agi("TT18"),agi("TT8"),agi("PER36"))
set.seed(35)
vioplot(na.exclude(contam_cors[contam_mat[['0']]]),na.exclude(contam_cors[contam_mat[['1']]]),na.exclude(contam_cors[contam_mat[['2:5']]]),na.exclude(contam_cors[contam_mat[['6:12']]]),
        add=F,names = c('0','1','2:5','6:12'),col = 'gray',rectCol = 'gray30',pchMed = 15,lwd = 1,drawRect = T,ylim = c(-1,1))
contam_mat[['0']]=vector()
contam_mat[['1']]=vector()
stripchart(lapply(contam_mat,function(x)contam_cors[x][contam_cors[x]<=0]), vertical = TRUE,group.names = names(contam_mat),
           method = "jitter",jitter = .2, add = TRUE, pch = 20,cex=.5, col = 'red')
abline(h=0,lty=2)
abline(v=2.5,lty=2)
points(x=rep(3,length(realmaternals[['2:5']])),y=contam_cors[realmaternals[['2:5']]],col='red',pch = dot)
text(x=rep(3,length(realmaternals[['2:5']])),y=contam_cors[realmaternals[['2:5']]],labels = c("FWA","FIS2"),col='black',pch = dot,pos = 4)
points(x=rep(4,length(realmaternals[['6:12']])),y=contam_cors[realmaternals[['6:12']]],col='red',pch = dot)
text(x=rep(4,length(realmaternals[['6:12']])),y=contam_cors[realmaternals[['6:12']]],labels = c("SDC","MRU1",'HDG9'),pos = 4,col='black')

points(x=rep(3,length(scmaternals[['2:5']])),y=contam_cors[scmaternals[['2:5']]],col='blue',pch = dot)
text(x=rep(3,length(scmaternals[['2:5']])),y=contam_cors[scmaternals[['2:5']]],labels = c("SHP1","SHP2","ELA1","TT1","TT5","TT16"),pos = 4,col='black')
points(x=rep(4,length(scmaternals[['6:12']])),y=contam_cors[scmaternals[['6:12']]],col='blue',pch = dot)
text(x=rep(4,length(scmaternals[['6:12']])),y=contam_cors[scmaternals[['6:12']]],labels = c("STK","BAN","TT10","TT3","TT19","TT4","TT7","TT18","TT8","PER36"),pos = 4,col='black')
cat('maternal_own_analysis:',length(which(contam_cors[names(which(allendomat>1))]>0)),
length(which(contam_cors[names(which(allendomat>1))]<=0))/length(which(allendomat>1)),
length(which(contam_cors[names(which(allendomat>1))]<=0)),'\n',sep='\n')

title(xlab = 'Number of Experiments Detected',ylab="Correlation with Seed Coat Contamination")
dev.off()
####################################################
contam_pat=list()
contam_pat[['0']]=names(contam_cors[!names(contam_cors)%in%names(allendopat)])
for(i in c('1','2:5','6:12')){
  contam_pat[[as.character(i)]]=names(allendopat[which(allendopat%in%eval(parse(text=i)))])
}
pdf("paternal_imprinting.pdf",7,7,useDingbats = F)
set.seed(35)
vioplot(na.exclude(contam_cors[contam_pat[['0']]]),na.exclude(contam_cors[contam_pat[['1']]]),na.exclude(contam_cors[contam_pat[['2:5']]]),na.exclude(contam_cors[contam_pat[['6:12']]]),
        add=F,names = c('0','1','2:5','6:12'),col = 'gray',rectCol = 'gray30',pchMed = 15,lwd = 1,drawRect = T,ylim = c(-1,1))
contam_pat[['0']]=vector()
contam_pat[['1']]=vector()
stripchart(lapply(contam_pat,function(x)contam_cors[x][contam_cors[x]<=0]), vertical = TRUE,group.names = names(contam_pat),
           method = "jitter",jitter = .2, add = TRUE, pch = 20,cex=.5, col = 'red')
abline(h=0,lty=2)
abline(v=2.5,lty=2)
points(x=c(4.08,4.07,4),y=contam_cors[c(agi("YUC10"),agi("HDG3"),agi("VIM5"))],col='red',pch = dot)
text(x=c(4.08,4.07,4),y=contam_cors[c(agi("YUC10"),agi("HDG3"),agi("VIM5"))],labels = c("YUC10","HDG3","VIM5"),pos = 4,col='black')
points(x=c(2.85,3),y=contam_cors[c(agi("NRPD1A"),agi("AGL92"))],col='red',pch = dot)
text(x=c(2.85,3),y=contam_cors[c(agi("NRPD1A"),agi("AGL92"))],labels = c("NRPD1A","AGL92"),col='black',pos = 4)
title(xlab = 'Number of Experiments Detected',ylab="Correlation with Seed Coat Contamination")
cat('paternal_own_analysis:',length(which(contam_cors[names(which(allendopat>1))]>0)),
    length(which(contam_cors[names(which(allendopat>1))]<0))/length(which(allendopat>1)),
    length(which(contam_cors[names(which(allendopat>1))]<=0)),'\n',sep='\n')
dev.off()

x=(na.exclude(contam_pvals[contam_pvals<.05]))
signifs=as.vector(x)
names(signifs)=as.character(names(x))
y=names(which(contam_cors>0))
z=names(which(contam_cors<0))
signif_positive=names(na.exclude(signifs[y]))
signif_negative=names(na.exclude(signifs[z]))


########################################################################################
# REPORTED IMPRINTED GENES
########################################################################################
reported_imprinting <- read.delim("../../contamination_paper/data_tables/reported_imprinting.txt", stringsAsFactors=FALSE)

rependomat=reported_imprinting$Total.Maternal[which(reported_imprinting$Total.Maternal>0)]
names(rependomat)=reported_imprinting$Gene.Identifier[which(reported_imprinting$Total.Maternal>0)]
rependopat=reported_imprinting$Total.Paternal[which(reported_imprinting$Total.Paternal>0)]
names(rependopat)=reported_imprinting$Gene.Identifier[which(reported_imprinting$Total.Paternal>0)]

dot=5
contam_mat=list()
contam_mat[['0']]=names(contam_cors[!names(contam_cors)%in%names(rependomat)])
for(i in c('1','2','3:5')){
  contam_mat[[as.character(i)]]=names(rependomat[which(rependomat%in%eval(parse(text=i)))])
}
pdf("reported_maternal_imprinting.pdf",7,7,useDingbats = F)
set.seed(32)
vioplot(na.exclude(contam_cors[contam_mat[['0']]]),na.exclude(contam_cors[contam_mat[['1']]]),na.exclude(contam_cors[contam_mat[['2']]]),na.exclude(contam_cors[contam_mat[['3:5']]]),
        add=F,names = c('0','1','2','3:5'),col = 'gray',rectCol = 'gray30',pchMed = 15,lwd = 1,drawRect = T,ylim = c(-1,1))
contam_mat[['0']]=vector()
contam_mat[['1']]=vector()
stripchart(lapply(contam_mat,function(x)contam_cors[x][contam_cors[x]<=0]), vertical = TRUE,group.names = names(contam_pat),
           method = "jitter",jitter = .2, add = TRUE, pch = 20,cex=.5, col = 'red')
stripchart(lapply(contam_mat,function(x)contam_cors[x][contam_cors[x]>0]), vertical = TRUE,group.names = names(contam_pat),
           method = "jitter",jitter = .2, add = TRUE, pch = 20,cex=.5, col = 'black')



abline(h=0,lty=2)

title(xlab = 'Number of Experiments Detected',main = "'Filtered' Maternally-Imprinted Genes Correlate with Seed Coat Contamination",ylab="Correlation with Severity of Contamination")
cat('maternal_reported:',length(which(contam_cors[names(which(rependomat>1))]>0)),
    length(which(contam_cors[names(which(rependomat>1))]<=0))/length(which(rependomat>1)),
    length(which(contam_cors[names(which(rependomat>1))]<=0)),'\n',sep='\n')
dev.off()
####################################################
contam_pat=list()
contam_pat[['0']]=names(contam_cors[!names(contam_cors)%in%names(rependopat)])
for(i in c('1','2','3:5')){
  contam_pat[[as.character(i)]]=names(rependopat[which(rependopat%in%eval(parse(text=i)))])
}
pdf("reported_paternal_imprinting.pdf",7,7,useDingbats = F)
plot(c(-10,-10),xlim=c(.5,4.5),ylim=c(-1,1),axes = F,xlab=NA,ylab=NA)
set.seed(35)
vioplot(na.exclude(contam_cors[contam_pat[['0']]]),na.exclude(contam_cors[contam_pat[['1']]]),na.exclude(contam_cors[contam_pat[['2']]]),na.exclude(contam_cors[contam_pat[['3:5']]]),
        add=F,names = c('0','1','2','3:5'),col = 'gray',rectCol = 'gray30',pchMed = 15,lwd = 1,drawRect = T,ylim = c(-1,1))

contam_pat[['0']]=vector()
contam_pat[['1']]=vector()
stripchart(lapply(contam_pat,function(x)contam_cors[x][contam_cors[x]<=0]), vertical = TRUE,group.names = names(contam_pat),
           method = "jitter",jitter = .2, add = TRUE, pch = 20,cex=.5, col = 'blue')
stripchart(lapply(contam_pat,function(x)contam_cors[x][contam_cors[x]>0]), vertical = TRUE,group.names = names(contam_pat),
           method = "jitter",jitter = .2, add = TRUE, pch = 20,cex=.5, col = 'black')
contam_pat[['0']]=names(contam_cors[!names(contam_cors)%in%names(rependopat)])
abline(h=0,lty=2)
title(xlab = 'Number of Experiments Detected',main = "Test for paternal Imprinting Correlates against Seed Coat Contamination",ylab="Correlation with Severity of Contamination")
cat('paternal_reported:',length(which(contam_cors[names(which(rependopat>1))]>0)),
    length(which(contam_cors[names(which(rependopat>1))]<=0))/length(which(rependopat>1)),
    length(which(contam_cors[names(which(rependopat>1))]<=0)),'\n',sep='\n')
dev.off()

############################################
# Correlation of tissue-enriched gene sets #
############################################

pdf('lc_gene_contamcor.pdf',6,4)

vioplot(na.exclude(contam_cors[!names(contam_cors)%in%unlist(enriched$lc)]),
        na.exclude(contam_cors[enriched$lc$EP]),
        na.exclude(contam_cors[enriched$lc$MCE]),
        na.exclude(contam_cors[enriched$lc$PEN]),
        na.exclude(contam_cors[enriched$lc$CZE]),
        na.exclude(contam_cors[enriched$lc$CZSC]),
        na.exclude(contam_cors[enriched$lc$GSC]),
        add=F,names = c('All genes',names(specifics$lc)[1:6]),col = 'gray',rectCol = 'gray30',pchMed = 15,lwd = 1,drawRect = T,ylim = c(-1,1))

abline(h=0,lty=2)
abline(v=c(2.5,5.5))
dev.off()

pdf('lc_whole_seed_bins.pdf',6,4)
x=rowMeans(belmonte_LCM[,84:85])
# x=x[names(x)[!names(x) %in% unique(unlist(lapply(enriched,unlist)))]]
x=x[grep(';',names(x),invert=T)]
x=sort(x)
avg_vals=list()
for(i in 1:5){
  avg_vals[[i]]=contam_cors[names(x)[(4172*(i-1)+1):(4172*i)]]
}

library(vioplot)
vioplot(na.exclude(avg_vals[[1]]),na.exclude(avg_vals[[2]]),na.exclude(avg_vals[[3]]),na.exclude(avg_vals[[4]]),na.exclude(avg_vals[[5]]),col = 'gray',ylim=c(-1,1))
wilcox.test(na.exclude(avg_vals[[1]]),na.exclude(avg_vals[[5]]))
abline(h=0,lty=2)
dev.off()




mean(contam_cors<=0,na.rm = T)
for(i in c('EP','MCE','PEN','CZE','CZSC','GSC')){
  cat(i,'stringent:',
      sum(contam_cors[enriched$lc[[i]]]<=0,na.rm = T),
      mean(contam_cors[enriched$lc[[i]]]<=0,na.rm = T),'\n',sep=' ')
  cat(i,'permissive:',
      sum(contam_cors[enriched$lc[[i]]]<=0 | contam_pvals[enriched$lc[[i]]]>=0.05,na.rm = T),
      mean(contam_cors[enriched$lc[[i]]]<=0 | contam_pvals[enriched$lc[[i]]]>=0.05,na.rm = T),'\n',sep=' ')
}
mean(contam_cors[enriched$lc$EP]<=0,na.rm = T)
mean(contam_cors[enriched$lc$MCE]<=0,na.rm = T)
mean(contam_cors[enriched$lc$PEN]<=0,na.rm = T)
mean(contam_cors[enriched$lc$CZE]<=0,na.rm = T)
mean(contam_cors[enriched$lc$CZSC]<=0,na.rm = T)
mean(contam_cors[enriched$lc$GSC]<=0,na.rm = T)

mean(contam_cors[enriched$lc$MCE]<=0 | contam_pvals[enriched$lc$MCE]>=0.05,na.rm = T)
mean(contam_cors[enriched$lc$PEN]<=0 | contam_pvals[enriched$lc$PEN]>=0.05,na.rm = T)
mean(contam_cors[enriched$lc$CZE]<=0 | contam_pvals[enriched$lc$CZE]>=0.05,na.rm = T)
mean(contam_cors[enriched$lc$CZSC]<=0 | contam_pvals[enriched$lc$CZSC]>=0.05,na.rm = T)
mean(contam_cors[enriched$lc$GSC]<=0 | contam_pvals[enriched$lc$GSC]>=0.05,na.rm = T)


for(i in c('EP','MCE','PEN','CZE','CZSC','GSC'))cat(i,ks.test(contam_cors,contam_cors[specifics$lc[[i]]],exact = T)$p.value,'\n',sep='\t')

cat('1_P_all',ks.test(contam_cors[!names(contam_cors)%in%names(allendopat)[which(allendopat%in%c(1))]],
                      contam_cors[names(allendopat)[which(allendopat%in%c(1))]],exact = T)$p.value,'\n',sep='\t')
cat('2-5_P_all',ks.test(contam_cors[!names(contam_cors)%in%names(allendopat)[which(allendopat%in%c(2:5))]],
                        contam_cors[names(allendopat)[which(allendopat%in%c(2:5))]],exact = T)$p.value,'\n',sep='\t')
cat('6-12_P_all',ks.test(contam_cors[!names(contam_cors)%in%names(allendopat)[which(allendopat%in%c(6:12))]],
                         contam_cors[names(allendopat)[which(allendopat%in%c(6:12))]],exact = T)$p.value,'\n',sep='\t')

cat('1_P_rep',ks.test(contam_cors[!names(contam_cors)%in%names(rependopat)[which(rependopat%in%c(1))]],
                      contam_cors[names(rependopat)[which(rependopat%in%c(1))]],exact = T)$p.value,'\n',sep='\t')
cat('2-5_P_rep',ks.test(contam_cors[!names(contam_cors)%in%names(rependopat)[which(rependopat%in%c(2))]],
                        contam_cors[names(rependopat)[which(rependopat%in%c(2))]],exact = T)$p.value,'\n',sep='\t')
cat('6-12_P_rep',ks.test(contam_cors[!names(contam_cors)%in%names(rependopat)[which(rependopat%in%c(3:5))]],
                         contam_cors[names(rependopat)[which(rependopat%in%c(3:5))]],exact = T)$p.value,'\n',sep='\t')

cat('1_M_all',ks.test(contam_cors[!names(contam_cors)%in%names(allendomat)[which(allendomat%in%c(1))]],
                      contam_cors[names(allendomat)[which(allendomat%in%c(1))]],exact = T)$p.value,'\n',sep='\t')
cat('1_M_all',ks.test(contam_cors[!names(contam_cors)%in%names(allendomat)[which(allendomat%in%c(2:5))]],
                      contam_cors[names(allendomat)[which(allendomat%in%c(2:5))]],exact = T)$p.value,'\n',sep='\t')
cat('1_M_all',ks.test(contam_cors[!names(contam_cors)%in%names(allendomat)[which(allendomat%in%c(6:12))]],
                      contam_cors[names(allendomat)[which(allendomat%in%c(6:12))]],exact = T)$p.value,'\n',sep='\t')

cat('1_M_rep',ks.test(contam_cors[!names(contam_cors)%in%names(rependomat)[which(rependomat%in%c(1))]],
                      contam_cors[names(rependomat)[which(rependomat%in%c(1))]],exact = T)$p.value,'\n',sep='\t')
cat('1_M_rep',ks.test(contam_cors[!names(contam_cors)%in%names(rependomat)[which(rependomat%in%c(2))]],
                      contam_cors[names(rependomat)[which(rependomat%in%c(2))]],exact = T)$p.value,'\n',sep='\t')
cat('1_M_rep',ks.test(contam_cors[!names(contam_cors)%in%names(rependomat)[which(rependomat%in%c(3:5))]],
                      contam_cors[names(rependomat)[which(rependomat%in%c(3:5))]],exact = T)$p.value,'\n',sep='\t')

