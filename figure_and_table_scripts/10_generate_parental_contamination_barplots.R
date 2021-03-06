splitRPM_all_p_values <- read.delim("enrichment_test_results/RPSM_p_values.txt",stringsAsFactors=F)
splitRPM_unique_p_values <- read.delim("enrichment_test_results/u_RPSM_p_values.txt",stringsAsFactors=F)


seedcoat_values=matrix(nrow=59,ncol=4)
colnames(seedcoat_values)=c('mat.all','mat.unique','pat.all','pat.unique')
rownames(seedcoat_values)=gsub('\\.pat','',grep('\\.pat',colnames(splitRPM_unique_p_values),value=T))
seedcoat_values[,1]=as.numeric(splitRPM_all_p_values[7,grep('\\.mat',colnames(splitRPM_all_p_values))])
seedcoat_values[,2]=as.numeric(splitRPM_unique_p_values[7,grep('\\.mat',colnames(splitRPM_unique_p_values))])
seedcoat_values[,3]=as.numeric(splitRPM_all_p_values[7,grep('\\.pat',colnames(splitRPM_all_p_values))])
seedcoat_values[,4]=as.numeric(splitRPM_unique_p_values[7,grep('\\.pat',colnames(splitRPM_unique_p_values))])

embryo_values=matrix(nrow=59,ncol=4)
colnames(embryo_values)=c('mat.all','mat.unique','pat.all','pat.unique')
rownames(embryo_values)=gsub('\\.pat','',grep('\\.pat',colnames(splitRPM_unique_p_values),value=T))
embryo_values[,1]=as.numeric(splitRPM_all_p_values[1,grep('\\.mat',colnames(splitRPM_all_p_values))])
embryo_values[,2]=as.numeric(splitRPM_unique_p_values[1,grep('\\.mat',colnames(splitRPM_unique_p_values))])
embryo_values[,3]=as.numeric(splitRPM_all_p_values[1,grep('\\.pat',colnames(splitRPM_all_p_values))])
embryo_values[,4]=as.numeric(splitRPM_unique_p_values[1,grep('\\.pat',colnames(splitRPM_unique_p_values))])

suspensor_values=matrix(nrow=59,ncol=4)
colnames(suspensor_values)=c('mat.all','mat.unique','pat.all','pat.unique')
rownames(suspensor_values)=gsub('\\.pat','',grep('\\.pat',colnames(splitRPM_unique_p_values),value=T))
suspensor_values[,1]=as.numeric(splitRPM_all_p_values[2,grep('\\.mat',colnames(splitRPM_all_p_values))])
suspensor_values[,2]=as.numeric(splitRPM_unique_p_values[2,grep('\\.mat',colnames(splitRPM_unique_p_values))])
suspensor_values[,3]=as.numeric(splitRPM_all_p_values[2,grep('\\.pat',colnames(splitRPM_all_p_values))])
suspensor_values[,4]=as.numeric(splitRPM_unique_p_values[2,grep('\\.pat',colnames(splitRPM_unique_p_values))])


endosperm_values=matrix(nrow=59,ncol=4)
colnames(endosperm_values)=c('mat.all','mat.unique','pat.all','pat.unique')
rownames(endosperm_values)=gsub('\\.pat','',grep('\\.pat',colnames(splitRPM_unique_p_values),value=T))
endosperm_values[,1]=as.numeric(splitRPM_all_p_values[4,grep('\\.mat',colnames(splitRPM_all_p_values))])
endosperm_values[,2]=as.numeric(splitRPM_unique_p_values[4,grep('\\.mat',colnames(splitRPM_unique_p_values))])
endosperm_values[,3]=as.numeric(splitRPM_all_p_values[4,grep('\\.pat',colnames(splitRPM_all_p_values))])
endosperm_values[,4]=as.numeric(splitRPM_unique_p_values[4,grep('\\.pat',colnames(splitRPM_unique_p_values))])



embryo_samples=c(1,2,3,4,5,6,8,13,14,15,17,19,20,21,22,23,27,36,40,44,46)
percentmat=as.numeric(gsub('%','',mapped_table[gsub('\\.','-',gsub('^u_','',rownames(endosperm_values)[embryo_samples])),"Percent.Maternal"]))
embryo_samples=embryo_samples[order(percentmat)]
endosperm_samples=c(9,11,28,12,7,25,26,16,18,37,38,41,42,45,47,48,49,50,51,52,53,54,55,56,57,58,59)
percentmat=as.numeric(gsub('%','',mapped_table[gsub('\\.','-',gsub('^u_','',rownames(endosperm_values)[endosperm_samples])),"Percent.Maternal"]))
endosperm_samples=endosperm_samples[order(percentmat)]
columns=c(1,3)

pdf('embryo_parental_contamination.pdf',8,8)
par(mar=c(4,5,2,2),mfrow=c(2,1))
barplot(t(-log10(embryo_values[embryo_samples,columns])),beside=T,las=2,col=c('firebrick2','deepskyblue2'),border=NA,ylim=c(0,30),axisnames=F)
abline(h=-log10(.05),lty=2)
abline(h=0,lty=1)
title(ylab='Embryo Proper Enrichment\n(-log10 p-value)')
legend('topleft',legend=c('Maternal','Paternal'),fill=c('firebrick2','deepskyblue2'),border=NA,box.col=NA,bg='white')
barplot(t(-log10(seedcoat_values[embryo_samples,columns])),beside=T,las=2,col=c('firebrick2','deepskyblue2'),border=NA,ylim=c(0,40))
abline(h=-log10(.05),lty=2)
abline(h=0,lty=1)
title(ylab='General Seed Coat Enrichment\n(-log10 p-value)')
legend('topleft',legend=c('Maternal','Paternal'),fill=c('firebrick2','deepskyblue2'),border=NA,box.col=NA,bg='white')
dev.off()

pdf('endosperm_parental_contamination.pdf',8,8)
par(mar=c(4,5,2,2),mfrow=c(2,1))
barplot(t(-log10(endosperm_values[endosperm_samples,columns])),beside=T,las=2,col=c('firebrick2','deepskyblue2'),border=NA,ylim=c(0,15),axisnames=F)
abline(h=-log10(.05),lty=2)
abline(h=0,lty=1)
title(ylab='Peripheral Endosperm Enrichment\n(-log10 p-value)')
abline(v=c(1,3),lty=3)
legend('topleft',legend=c('Maternal','Paternal'),fill=c('firebrick2','deepskyblue2'),border=NA,box.col=NA)
barplot(t(-log10(seedcoat_values[endosperm_samples,columns])),beside=T,las=2,col=c('firebrick2','deepskyblue2'),border=NA,ylim=c(0,50))
abline(h=-log10(.05),lty=2)
abline(h=0,lty=1)
title(ylab='General Seed Coat Enrichment\n(-log10 p-value)')
abline(v=c(1,3),lty=3)
legend('topleft',legend=c('Maternal','Paternal'),fill=c('firebrick2','deepskyblue2'),border=NA,box.col=NA)
dev.off()