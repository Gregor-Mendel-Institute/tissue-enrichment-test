sourcedir="C:/Users/schon.admin/Desktop/TPC_revisions/datasets_generated/"
setwd(sourcedir)
snps_per_gene=read.table(paste(sourcedir,'/snps_per_gene.tsv',sep=''),row.names = 1,header = T,stringsAsFactors = F)


# for(prefix in c('u_','4MM_','')){
for(prefix in c('')){
  genecounts=c('24cell_KypCol_11','24cell_kypCol_12','24cell_LerCol_11','24cell_LerCol_7','24cell_SC_2','2cell_ColxCvi','2cell_CvixCol','8cell_ColxCvi','8cell_CvixCol','8cell_pilot_ColxCvi','8cell_pilot_CvixCol','Col-dmexLer-WT_Endosperm','Col-LCMxLer-LCM_Embryo','Col-LCMxLer-LCM_Endosperm','Col-WTxLer-WT_Embryo','Col-WTxLer-WT_Endosperm','Col-WTxLer-WT_full-Endosperm','ColxCol_embryo','ColxCol_endosperm','ColxCvi_embryo','ColxCvi_endosperm_1','ColxCvi_endosperm_2','ColxCvi_endosperm_3','ColxCvi_whole_seed','ColxLer_endosperm_2','ColxLer_endosperm_4','ColxLer_endosperm_5','CvixCol_embryo','CvixCol_endosperm_1','CvixCol_endosperm_2','CvixCol_endosperm_3','CvixCol_whole_seed','CvixCvi_embryo','CvixCvi_endosperm','CvixLer_embryo','CvixLer_endosperm','CvixLer_endosperm_2','CvixLer_endosperm_3','dcl1_globular','eg_ColxCvi','eg_CvixCol','gehring_ColLer_emb','gehring_ColLer_endo','gehring_LerCol_emb','gehring_LerCol_endo','gl_LerCol_4','Ler-fiexCol-WT_Endosperm','Ler-WTxCol-met1_Endosperm','Ler-WTxCol-WT_Embryo','Ler-WTxCol-WT_Endosperm','Ler-WTxCol-WT_full-Endosperm','LerxCol_endosperm_2','LerxCol_endosperm_4','LerxCol_endosperm_5','LerxCvi_embryo','LerxCvi_endosperm','LerxCvi_endosperm_2','LerxCvi_endosperm_3','LerxLer_embryo','LerxLer_endosperm','wt_globular')
  genecounts = paste(prefix,genecounts,sep='')
  guide_table <- read.delim(paste(sourcedir,"sample_data_table.txt",sep=''), stringsAsFactors=FALSE)
  rownames(guide_table)=guide_table[,1]
  rownames(guide_table)=paste(prefix,rownames(guide_table),sep='')
  genecounts=genecounts[genecounts%in%rownames(guide_table)]
  cross_table=list()
  for(i in genecounts){
    if(i==paste(prefix,'dcl1_globular',sep='')){
      countsfile = grep(paste(gsub(prefix,'','GaD_dcl1'),'.',sep=''),list.files(paste(sourcedir,prefix,"genecounts/",sep="")),value = T,fixed = T)
    }else{
      if(i==paste(prefix,'wt_globular',sep='')){
        countsfile = grep(paste(gsub(prefix,'','GaD_wt'),'.',sep=''),list.files(paste(sourcedir,prefix,"genecounts/",sep="")),value = T,fixed = T)
      }else{
        countsfile = grep(paste(gsub(prefix,'',i),'.',sep=''),list.files(paste(sourcedir,prefix,"genecounts/",sep="")),value = T,fixed = T)
      }
    }
    tmp=read.delim(paste(c(sourcedir,prefix,"genecounts/",countsfile),collapse=""),header = F,sep="\t",row.names = 1)
    tmpnames=sort(unique(gsub("(.+)(col|cvi|ler)","\\1",rownames(tmp))))
    parents=names(sort(table(substring(rownames(tmp),10,13)),decreasing = T)[1:2])

    cross_table[[i]]=matrix(nrow=length(tmpnames),ncol=2)
    rownames(cross_table[[i]])=tmpnames
    colnames(cross_table[[i]])=parents
    cross_table[[i]][,parents[1]]=tmp[paste(tmpnames,parents[1],sep=""),1]
    cross_table[[i]][,parents[2]]=tmp[paste(tmpnames,parents[2],sep=""),1]
  }
  cat(prefix,'')
  cat(c('Sample','mat.reads','pat.reads','percent.mat','mat.genes','pat.genes'),sep='\t')
  cat('\n')
  all_counts=matrix(nrow=nrow(cross_table[[1]]))
  for(i in genecounts){
    allreads=(colSums(cross_table[[i]],na.rm=T))
    maternal=strsplit(tolower(guide_table[i,'Maternal.Ecotype']),'-')[[1]][1]
    paternal=strsplit(tolower(guide_table[i,'Paternal.Ecotype']),'-')[[1]][1]
    parents=names(allreads)
    if(maternal==paternal){
      paternal=parents[which(parents!=maternal)]
    }
    cat(c(gsub('^(u|4MM)_','',gsub('\\.counts','',gsub('(SRX|GSM)[0-9]+_','',i))),
          allreads[maternal],
          allreads[paternal],
          round(allreads[maternal]/sum(allreads),4)*100,
          length(which(cross_table[[i]][,maternal]>0)),
          length(which(cross_table[[i]][,paternal]>0))),sep='\t')
    cat('\n')
    tablenames = colnames(all_counts)
    if(is.null(tablenames)){
      tablenames=1
    }
    all_counts = cbind(all_counts, cross_table[[i]][,maternal],cross_table[[i]][,paternal])
    colnames(all_counts) = c(tablenames,paste(i,c('mat','pat'),sep='.'))
  }
  all_counts=all_counts[,2:ncol(all_counts)]
  all_counts=all_counts[rownames(snps_per_gene),]
  RPSM=all_counts
  for(i in 1:ncol(RPSM)){
    samplename=strsplit(colnames(RPSM)[i],'.',fixed = T)[[1]][1]
    ecotypes=sort(c(substring(guide_table[samplename,"Maternal.Ecotype"],1,3),substring(guide_table[samplename,"Paternal.Ecotype"],1,3)))
    if(ecotypes[1]==ecotypes[2]){
      if(ecotypes[1]=='Col'){
        ecotypes[2]='Ler'
      }else{
        if(ecotypes[1]=='Ler'){
          ecotypes[1]='Col'
        }else{
          if(ecotypes[1]=='Cvi'){
            ecotypes[1]='Col'
          }
        }
      }
    }
    eco = paste(ecotypes,collapse='_')
    cat(c(samplename,eco),'\n',sep='\t')
    snps=snps_per_gene[,eco]
    RPSM[,i]=RPSM[,i]/snps
    RPSM[is.nan(RPSM[,i]) | is.infinite(RPSM[,i]),i] = NA
  }
  for(i in seq(1,ncol(RPSM),2)){
    column1=RPSM[,i]/sum(colSums(RPSM[,i:(i+1)],na.rm = T))*10^6
    column2=RPSM[,i+1]/sum(colSums(RPSM[,i:(i+1)],na.rm = T))*10^6
    RPSM[,c(i,i+1)]=cbind(column1,column2)
  }
  # write.table(all_counts,paste(prefix,'all_counts.tsv',sep=''),sep='\t',quote=F)
  # write.table(RPSM,paste(prefix,'RPSM.tsv',sep=''),sep='\t',quote=F)
}


barplot(colSums(RPSM,na.rm = T),space=rep(c(1,0),ncol(RPSM)/2),names.arg = NA,col=c('firebrick2','deepskyblue2'))
axis(1,seq(2,by = 3,length.out = ncol(RPSM)/2),labels = unique(gsub('\\.[mp]at','',colnames(RPSM))),las=2)


rm(tmpnames,tmp,parents,x,uniquereads,allreads)
