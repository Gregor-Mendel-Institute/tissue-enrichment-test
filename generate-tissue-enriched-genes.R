##########################
# STATISTICAL THRESHOLDS #
#############################################################################################################################

minimum_fold <- 8              # Smallest fold change to be considered 'enriched' in a tissue
testing_correction <- 'fdr'    # P-value adjustment for multiple testing correction
maximum_pvalue <- 0.001        # Largest adjusted P-value to be considered significant in ANOVA
flanking <- TRUE               # If TRUE, buffers temporal changes by including adjacent developmental timepoints in analysis
disjoint_sets <- FALSE          # If TRUE, ensures that all tissue-enriched gene sets for each stage are disjoint (non-overlapping)
reference_name <- 'reference_atlas'

#############################################################################################################################
# FUNCTIONS #
#############################################################################################################################
# Takes a list of sets and returns the same list with all elements shared between any sets removed
disjoint <- function(input_list){ 
  disjoint_list <- input_list
  for(i in names(input_list)){
    for(j in names(input_list)[names(input_list)!=i]){
      disjoint_list[[i]] <- disjoint_list[[i]][!disjoint_list[[i]]%in%disjoint_list[[j]]]
    }
  }
  return(disjoint_list)
}

#############################################################################################################################
# LOADING ENVIRONMENT #
#############################################################################################################################
args = commandArgs(trailingOnly=TRUE)
if(length(args)==1){
	reference_name = gsub('^datasets/','',gsub('\\.(tsv|description)$','',args[1]))
}
if(grepl('/',reference_name)){
	stop("Reference table must be located in the local directory 'datasets'. Please check setup instructions (See README for details)")
}

cat('Loading reference atlas...\n')
reference_atlas <- read.delim(paste("datasets/",reference_name,".tsv",sep=''),header = T,stringsAsFactors = F,comment.char = '#')
tissues <- scan("tissues.txt",'character',sep='\n',comment.char = '#')
tissues <- strsplit(tissues,'\t')
names(tissues) <- unlist(lapply(tissues,function(x)x[1]))
tissues <- unlist(lapply(tissues,function(x)x[2]))
tissue_cols <- gsub('^(.+?)_?(.+)\\..*$','\\1',colnames(reference_atlas))
if(!all(tissue_cols%in%names(tissues))){
  cat("WARNING: Unrecognized tissue names! Ignoring samples:\n\t")
  ignr <- colnames(reference_atlas)[!tissue_cols%in%names(tissues)]
  cat(ignr,sep='\t')
  reference_atlas <- reference_atlas[,!colnames(reference_atlas)%in%ignr, drop=FALSE]
}
excluded_tissues <- names(tissues)[grep('//exclude',tissues)]
if(length(excluded_tissues)>0){
  for(i in 1:length(excluded_tissues)){
    tissue_cols <- gsub('^(.+?)_?(.+)\\..*$','\\1',colnames(reference_atlas))
    tissues_to_exclude <- grep(excluded_tissues[i],tissue_cols)
    cat(paste('Removed',length(tissues_to_exclude),gsub('//exclude','',tissues[excluded_tissues[i]]),'samples from reference atlas.\n',sep=' '))
    x <- 1:ncol(reference_atlas)
    reference_atlas <- reference_atlas[,x[!x%in%tissues_to_exclude], drop=FALSE]
  }
}
timepoints <- scan("timepoints.txt",'character',sep='\n',comment.char = '#')
timepoints <- strsplit(timepoints,'\t')
names(timepoints) <- unlist(lapply(timepoints,function(x)x[1]))
timepoints <- unlist(lapply(timepoints,function(x)x[2]))
timepoint_cols <- gsub('^(.+?)_?(.+)\\..*$','\\2',colnames(reference_atlas))
if(!all(timepoint_cols%in%names(timepoints))){
  cat("WARNING: Unrecognized timepoint names! Ignoring samples:\n\t")
  ignr <- colnames(reference_atlas)[!timepoint_cols%in%names(timepoints)]
  cat(ignr,sep='\t')
  reference_atlas <- reference_atlas[,!colnames(reference_atlas)%in%ignr, drop=FALSE]
}
excluded_timepoints <- names(timepoints)[grep('//exclude',timepoints)]
if(length(excluded_timepoints)>0){
  for(i in 1:length(excluded_timepoints)){
    timepoint_cols <- gsub('^(.+?)_?(.+)\\..*$','\\2',colnames(reference_atlas))
    timepoints_to_exclude <- grep(excluded_timepoints[i],timepoint_cols)
    cat(paste('Removed',length(timepoints_to_exclude),gsub('//exclude','',timepoints[excluded_timepoints[i]]),'samples from reference atlas.\n',sep=' '))
    x=1:ncol(reference_atlas)
    reference_atlas <- reference_atlas[,x[!x%in%timepoints_to_exclude], drop=FALSE]
  }
}
reference_atlas <- reference_atlas[grep(';',rownames(reference_atlas),invert=TRUE),]
cat("Reference atlas prepared.\nNumber of samples detected at each stage:\n")

#############################################################################################################################
# PERFORMING GENE SELECTION #
#############################################################################################################################

tissue_cols <- gsub('^(.+?)_?(.+)\\..*$','\\1',colnames(reference_atlas))
timepoint_cols <- gsub('^(.+?)_?(.+)\\..*$','\\2',colnames(reference_atlas))
rownames(reference_atlas) <- toupper(rownames(reference_atlas))
for(i in names(timepoints)){
  cat(timepoints[i],'\n',sep='')
  tbl <- table(tissue_cols[grep(i,timepoint_cols)])
  if(length(tbl)==0){
    cat('\tNo samples')
  }else{
    cat(paste('\t',tissues[sort(names(tbl))],': ',tbl[sort(names(tbl))],'\n',sep=''))
  }
}

cat('\nIdentifying tissue-enriched gene expression...\n')
enriched <- list()
minimum_fold <- log2(minimum_fold)
for(stage in 1:length(timepoints)){
  enriched[[names(timepoints)[stage]]]=list()
  if(flanking){
    samples <- which(timepoint_cols%in%as.character((stage-1):(stage+1)))
  }else{
    samples <- which(timepoint_cols==as.character(stage))
  }
  subreference_atlas <- reference_atlas[,samples, drop=FALSE]
  for(tissue in names(tissues)){
    tissuereps <- grep(paste('^',tissue,'_',sep=''),colnames(subreference_atlas),ignore.case = F,value = T)
    nottissuereps <- colnames(subreference_atlas)[!colnames(subreference_atlas)%in%tissuereps]
    if(length(tissuereps)==0){
      next
    }
    cat(paste(tissuereps,collapse <- ','),'\nvs.\n',paste(nottissuereps,collapse=','),'\n',sep='')
    foldchanges <- log2(rowMeans(subreference_atlas[,tissuereps, drop=FALSE]+.1)/rowMeans(subreference_atlas[,nottissuereps, drop=FALSE]+.1))
    spec_fc <- names(which(foldchanges>=minimum_fold))
    subsubreference_atlas <- subreference_atlas[spec_fc,c(tissuereps,nottissuereps)]
    pvals <- apply(subsubreference_atlas,1,function(x){
      dt <- data.frame(values=as.numeric(x),group=c(rep('in',length(tissuereps)),rep('out',length(nottissuereps))),stringsAsFactors = F)
      dt.aov <- aov(dt$values ~ dt$group)
      return(as.numeric(unlist(summary(dt.aov))['Pr(>F)1']))
    })
    pvals <- p.adjust(pvals,method = testing_correction)
    pvalhits <- names(pvals[pvals<maximum_pvalue])
    cat(length(pvalhits),'\n\n')
	if(disjoint_sets){
	pvalhits <- disjoint(pvalhits)
	}
    enriched[[names(timepoints)[stage]]][[tissue]]<-pvalhits
  }
}

cat('',file = 'all_tested_genes.txt',append = F)
cat(paste(rownames(reference_atlas),collapse='\n'),file = 'all_tested_genes.txt',append=T)
cat("All genes in the reference atlas listed in 'all_tested_genes.txt'.\n")

cat('',file = 'enriched_genes.txt',append = F)
for(i in names(timepoints)){
  for(j in names(enriched[[i]])){
    cat(paste(paste(i,j,sep='_'),paste(enriched[[i]][[j]],collapse=','),sep='\t'),'\n',sep='',file = 'enriched_genes.txt',append=T)
  }
}
cat("Tissue-enriched transcripts saved in 'enriched_genes.txt'.")