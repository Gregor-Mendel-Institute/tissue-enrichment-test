setwd('C:/Users/schon.admin/Desktop/TPC_revisions/datasets_generated/')
enrichment_scores=matrix(nrow=7,ncol=1)
rownames(enrichment_scores)=c('EP','SUS','MCE','PEN','CZE','CZSC','GSC')
for(i in c('array_master','day_2007','reference_atlas','rsem_TPM','xiang_tableS1','RPSM')){
  tmp_table=read.table(paste('enrichment_test_results/',i,'_p_values.txt',sep=''),sep='\t',header=T,stringsAsFactors = F)
  enrichment_scores=cbind(enrichment_scores,-log10(tmp_table))
}
enrichment_scores=enrichment_scores[,2:ncol(enrichment_scores)]
write.table(t(enrichment_scores),'enrichment_scores.tsv',quote=F,sep='\t')

data_table=read.delim(paste(sourcedir,'all_samples_in_study.tsv',sep=''),stringsAsFactors = F,comment.char = '#',sep='\t')
rownames(data_table)=data_table[,1]

mapped_table=read.delim(paste(sourcedir,'read_mapping_statistics.tsv',sep=''),stringsAsFactors = F,comment.char = '#',sep='\t')
rownames(mapped_table)=mapped_table[,1]

######### PERCENT #########
pdf('scatter_embryo.pdf',6,6,useDingbats = F)
x=rownames(data_table)[which(data_table$Tissue.Type=='embryo'
        & data_table$Transcriptome.Type=='mRNA-seq'
        #& data_table$Maternal.Genotype =='WT'
        #& data_table$Paternal.Genotype =='WT'
        # & data_table$Corresponding.Developmental.Stage%in%c('preglobular','globular')
        & data_table$Maternal.Ecotype != data_table$Paternal.Ecotype)]
percentmat=as.numeric(gsub('%','',mapped_table[x,"Percent.Maternal"]))
colors=rep('black',length(x))
colors[which(data_table[x,"General.Seed.Coat"]>1.3)]='red'
colors[which(rowSums(data_table[x,c(9,11)]=='WT')<2)]='gray'

xvalues=data_table[x,"General.Seed.Coat"]
# xvalues=log2((rowSums(data_table[x,18:19])+.5)/(rowSums(data_table[x,13:14],na.rm = T)+.5))

xlim=c(0,50)
ylim=c(40,100)

plot(xvalues,percentmat,pch=20,
     xlab='-log10 p-value of Distal Seed Coat Enrichment',ylab='Percent of Maternally-Derived Reads',
     main='Effect of Seed Coat Contamination on Transcriptome Maternal Dominance',cex.lab=.8,cex.main=.9,
     xlim=xlim,col=colors,
     ylim=ylim,las=1)
text(xvalues,percentmat,labels = data_table[x,1],cex=.5,pos = 4)
xvalues=xvalues[which(rowSums(data_table[x,c(9,11)]=='WT')==2)]
percentmat=percentmat[which(rowSums(data_table[x,c(9,11)]=='WT')==2)]
y=cor.test(xvalues,percentmat,exact = T)
text(xlim[1],ylim[2],paste("Pearson's Correlation: ",signif(y$estimate,digits=3),sep=''),pos = 4)
text(xlim[1],ylim[2]-5,paste("p-value: ",format(signif(y$p.value,digits=3),scientific = T),sep=''),pos = 4)
abline(lm(percentmat~xvalues))
abline(h=50,lty=2)
dev.off()

pdf('scatter_endosperm.pdf',6,6,useDingbats = F)
x=which(data_table$Tissue.Type%in%c('endosperm','whole seed')
        & data_table$Transcriptome.Type=='mRNA-seq'
        #& data_table$Maternal.Genotype =='WT'
        #& data_table$Paternal.Genotype =='WT'
        & data_table$Maternal.Ecotype != data_table$Paternal.Ecotype)
colors=rep('black',length(x))
colors[which(data_table[x,"General.Seed.Coat"]>1.3)]='red'
colors[which(rowSums(data_table[x,c(9,11)]=='WT')<2)]='gray'
colors[which(x==grep('CvixLer_endosperm_3',data_table[,1]))]='gray'

xvalues=data_table$General.Seed.Coat[x]
# xvalues=log2((rowSums(data_table[x,18:19])+.5)/(rowSums(data_table[x,15:17])+.5))

percentmat=as.numeric(gsub('%','',mapped_table[rownames(data_table)[x],"Percent.Maternal"]))
xlim=c(0,50)
ylim=c(40,100)
plot(xvalues,percentmat,
     xlab='-log10 p-value of Distal Seed Coat Enrichment',ylab='Percent of Maternally-Derived Reads',
     main='Effect of Seed Coat Contamination on Transcriptome Maternal Dominance',pch=20,cex.lab=.8,cex.main=.9,
     xlim=xlim,col=colors,
     ylim=ylim,las=1)
text(xvalues,percentmat,labels = data_table[x,1],cex=.5,pos = 4)

x=which(data_table$Tissue.Type%in%c('endosperm','whole seed')
        & data_table$Transcriptome.Type=='mRNA-seq'
        & data_table$Maternal.Genotype =='WT'
        & data_table$Paternal.Genotype =='WT'
        & data_table$Maternal.Ecotype != data_table$Paternal.Ecotype)
x=x[!x%in%grep('CvixLer_endosperm_3',data_table[,1])]
xvalues=data_table$General.Seed.Coat[x]
percentmat=as.numeric(gsub('%','',mapped_table[rownames(data_table)[x],"Percent.Maternal"]))
# xvalues=log2((rowSums(data_table[x,18:19])+.5)/(rowSums(data_table[x,15:17])+.5))

y=cor.test(xvalues,percentmat,exact = T)
text(xlim[1],ylim[2],paste("Pearson's Correlation: ",signif(y$estimate,digits=3),sep=''),pos = 4)
text(xlim[1],ylim[2]-5,paste("p-value: ",format(signif(y$p.value,digits=3),scientific = T),sep=''),pos = 4)
abline(lm(percentmat~xvalues))
abline(h=2/3*100,lty=2)

dev.off()

