reported_but_greater_than_zero = sort(contam_cors[names(which(contam_cors[names(which(rependomat>1))]>0))],decreasing = T)
reported_and_not_cor = names(which(contam_cors[names(which(rependomat>1))]>0))


(reported_but_greater_than_zero)


xin2013_mat=scan('xin_2013_mat.txt','character',sep='\n')
xin2013_pat=scan('xin_2013_pat.txt','character',sep='\n')

xin2013_mat=sort(unique(xin2013_mat))
xin2013_pat=sort(unique(xin2013_pat))
vioplot(contam_cors[xin2013_mat],contam_cors[xin2013_pat])
abline(h=0)
intersect(names(which(contam_cors[xin2013_pat]<=0)),
          names(which(contam_cors[names(which(allendopat>1))]<=0)))

TAIRnames(intersect(names(allendopat),(xin2013_pat)))

chosen_samples=data_table[data_table[,"Tissue.Type"]%in%c('endosperm','whole seed') & data_table[,"Maternal.Ecotype"]!=data_table[,"Paternal.Ecotype"],1]
RPSM_endosperm=RPSM[,gsub('^X','',gsub('\\.','-',gsub('\\.[mp]at','',colnames(RPSM))))%in%chosen_samples]
colnames(RPSM_endosperm)=gsub('^X','',gsub('\\.','-',colnames(RPSM_endosperm)))
nms=names(sort(contam_ratio[names(contam_ratio)%in%chosen_samples]))
RPSM_endosperm=RPSM_endosperm[,paste(rep(nms,each=2),c('mat','pat'),sep='-')]
RPSM_endosperm=RPSM_endosperm[,grep('(CvixLer|LerxCvi)_endosperm_3',colnames(RPSM_endosperm),invert = T)]

endosperm_barplot=function(gene){
  # par(mfrow=c(2,1))
  par(mar=c(16,4,2,2))
  colors=rep(c('firebrick2','deepskyblue2'),ncol(RPSM_endosperm)/2)
  colors[c(grep('ColxCol',colnames(RPSM_endosperm)),grep('CvixCvi',colnames(RPSM_endosperm)),grep('LerxLer',colnames(RPSM_endosperm)))]='gray'
  barplot(as.numeric(RPSM_endosperm[gene,]),col=colors,space=rep(c(1,0),ncol(RPSM_endosperm)/2),las=2)
  axis(1,seq(2,by = 3,length.out = ncol(RPSM_endosperm)/2),labels = unique(gsub('-[mp]at','',colnames(RPSM_endosperm))),las=2)
  print(t.test(RPSM_endosperm[gene,seq(1,ncol(RPSM_endosperm),2)]/2,RPSM_endosperm[gene,seq(2,ncol(RPSM_endosperm),2)],paired = T))
}

barplot(colMeans(RPSM_endosperm,na.rm = T),space=rep(c(1,0),ncol(RPSM_endosperm)/2),names.arg = NA)
axis(1,seq(2,by = 3,length.out = ncol(RPSM_endosperm)/2),labels = unique(gsub('-[mp]at','',colnames(RPSM_endosperm))),las=2)
