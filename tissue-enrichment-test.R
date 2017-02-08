################
# USER OPTIONS #
#############################################################################################################################
filename <- 'nodine_2012_embryos' # Name of the .tsv file containing transcriptomes to be tested. Try also: 'reference_atlas'
make_heatmaps <- TRUE             # If TRUE, generates a heatmap of p-values for all samples. Requires package 'pheatmap'
row_breaks <- c(2,5)              # A vector containing row numbers where a gap should be placed on the heatmap
col_breaks <- FALSE               # A vector containing column numbers where a gap should be placed on the heatmap
tissue_colors=c('#609A33','#C0D84E','#C56DAC','#F9C9DE','#F26D7E','#2A89B6','#66C5E9') # For CDF plots; one color per tissue
names(tissue_colors)=c("EP","SUS","MCE","PEN","CZE","CZSC","GSC")                      # All tissues used in the test

#############################################################################################################################
# FUNCTIONS #
#############################################################################################################################

trim <- function(x){
  gsub("^\\s+|\\s+$", "", x)
}

test_all_tissues <- function(seqdata,test_timepoint,colorset=tissue_colors){
# Performs a two-sided Wilcoxon rank-sum  on a population of genes and their expression values for each tissue-enriched gene set.
# Required inputs: 
#          seqdata- a named numeric vector containing an expression value for each gene
#          test_timepoint- developmental timepoint to compare to. Options in 'timepoints.txt'
	wilcox.pvals=list()
	percentile_seqdata <- sapply(seqdata,function(x)sum(seqdata<=x)/length(seqdata)*100)
    plot(c(-10,-10),xlim=c(0,100),ylim=c(0,1),las=1,xlab='Percentile',ylab='Cumulative Frequency',main='')
    for(i in 1:length(enriched[[test_timepoint]])){
      tissue <- names(enriched[[test_timepoint]])[i]
      E <- enriched[[test_timepoint]][[i]]
      not_E <- all_tested_genes[all_tested_genes %in% names(percentile_seqdata) & !all_tested_genes %in% E]
      E <- E[E%in%names(percentile_seqdata)]
      not_E <- not_E[not_E%in%names(percentile_seqdata)]
      values_E <- as.numeric(percentile_seqdata[E])
      values_not_E <- as.numeric(percentile_seqdata[not_E])
      if(mean(values_E) <= mean(values_not_E)){
        test_pvalue <- 1
      }else{
        test_pvalue <- wilcox.test(x = values_E,y = values_not_E,alternative = 'two.sided',exact = F)$p.value
      }
	  if(test_pvalue < 10^-50){
		test_pvalue <- 10^-50
	  }
      wilcox.pvals[[tissue]]=test_pvalue
      lines(ecdf(values_E),verticals=T,do.points=F,lwd=2,col=tissue_colors[tissue],cex=.5)
    }
	lines(ecdf(percentile_seqdata),verticals=T,do.points=F,lwd=2)
	output <- as.numeric(wilcox.pvals)
	names(output) <- names(wilcox.pvals)
	return(output)
}

perform_enrichment_test <- function(data_table,columns,timepoints,gaps_row,gaps_col,output_name,make_heatmaps){
  outputdir <- 'results'
  pvalue_matrix <- matrix(nrow=length(unique(unlist(lapply(enriched,names)))),ncol=length(columns))
  rownames(pvalue_matrix) <- unique(unlist(lapply(enriched,names)))
  for(i in 1:length(columns)){
    cat("Performing tissue-enrichment test on ",columns[i],'\n')
    pdf(paste(outputdir,'/',columns[i],".pdf",sep=""),4,6)
    seqdata <- data_table[,columns[i]]
    seqdata <- as.numeric(trim(seqdata))
    names(seqdata) <- rownames(data_table)
    seqdata <- seqdata[!(is.na(seqdata))]
	seqdata <- seqdata[names(seqdata) %in% all_tested_genes]
    pvals <- test_all_tissues(seqdata,test_timepoint = timepoints[i])
    pvalue_matrix[names(pvals),i] <- as.numeric(pvals)
	text(x=-5,y=.975,labels='Enrichment Score',cex=.8,pos=4)
	legend(x=-5,y=.975,legend=paste(names(pvals),round(-log10(pvals),1),sep=': '),fill=tissue_colors[names(pvals)],bty='n',cex=.8,border=NA)
    title(main=columns[i],cex=.8)
    dev.off()
  }
  colnames(pvalue_matrix) <- columns
  cat('Storing p-values in ',paste(outputdir,paste(output_name,'p_values.txt',sep="_"),sep='/'),'\n',sep='')
  write.table(pvalue_matrix,paste(outputdir,paste(output_name,'p_values.txt',sep="_"),sep='/'),sep="\t",quote=F)
  if(make_heatmaps){
    cat('Saving heatmap to ',paste(outputdir,paste(output_name,"_heatmap.pdf",sep=''),sep='/'),'\n',sep='')
    breaksList = c(0,-log10(0.05),-log10(0.01),3,4,5,10,15,20,30,40,50)
    pdf(paste(outputdir,paste(output_name,"_heatmap.pdf",sep=''),sep='/'),7,4,onefile = F)
    pheatmap(-log10(pvalue_matrix[names(tissues),]),cluster_rows = F,cluster_cols = F,scale = 'none', breaks = breaksList,color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)),
             border_color = 'black',gaps_col = gaps_col,
             gaps_row = gaps_row,main="Tissue-Specific Transcript Enrichment",
             labels_row = tissues,
             labels_col = colnames(pvalue_matrix),fontsize = 6)
    dev.off()
  }
}

#############################################################################################################################
# LOADING ENVIRONMENT #
#############################################################################################################################
args = commandArgs(trailingOnly=TRUE)
if(length(args)==1){
	filename = gsub('^datasets/','',gsub('\\.(tsv|description)$','',args[1]))
}
if(grepl('/',filename)){
	stop("Data table must be located in the local directory 'datasets'. Please make a .tsv file and .description file of the dataset. (See README for details)")
}

if(make_heatmaps==TRUE){
	for(package in c('RColorBrewer','pheatmap')){
		if(!require(package,lib.loc = 'Rpackages',character.only = TRUE)){
			install.packages(package,lib = 'Rpackages',repos = "http://cran.us.r-project.org",dependencies = NA)
			library(package,lib.loc = 'Rpackages',character.only = TRUE)
		}
	}
}

tissues_file <- scan("tissues.txt",'character',sep='\n',comment.char = '#')
tissues_file <- strsplit(tissues_file,'\t')
tissues <- unlist(lapply(tissues_file,function(x)x[2]))
names(tissues) <- unlist(lapply(tissues_file,function(x)x[1]))
tissues <- tissues[grep('//exclude',tissues,invert=T)]

timepoints <- scan("timepoints.txt",'character',sep='\n',comment.char = '#')
timepoints <- strsplit(timepoints,'\t')
names(timepoints) <- unlist(lapply(timepoints,function(x)x[1]))
timepoints <- unlist(lapply(timepoints,function(x)x[2]))

enriched_genes <- scan("enriched_genes.txt",'character',sep='\n')
enriched_samples <- unlist(lapply(strsplit(enriched_genes,'\t'),function(x)x[1]))
enriched_genelist <- unlist(lapply(strsplit(enriched_genes,'\t'),function(x)x[2]))
enriched <- list()
for(i in 1:length(enriched_genes)){
  timepoint <- gsub('^([0-9]+)_(.+)$','\\1',enriched_samples[i])
  tissue <- gsub('^([0-9]+)_(.+)$','\\2',enriched_samples[i])
  enriched[[timepoint]][[tissue]] <- unlist(strsplit(enriched_genelist[i],','))
}
rm(enriched_genes,enriched_samples,enriched_genelist)

transcriptomes <- read.delim(paste('datasets/',filename,'.tsv',sep=''),header = T,stringsAsFactors = F,row.names=1)
rownames(transcriptomes) <- toupper(trim(rownames(transcriptomes)))
colnames(transcriptomes) <- trim(colnames(transcriptomes))

transcriptome_descriptions <- read.delim(paste('datasets/',filename,'.description',sep=''),header = T,stringsAsFactors = F)
rownames(transcriptome_descriptions) <- transcriptome_descriptions[,1]
samples <- trim(rownames(transcriptome_descriptions))

samples <- gsub('-','\\.',samples)
appended_x <- gsub('^X','',colnames(transcriptomes))
colnames(transcriptomes)[appended_x %in% samples] <- appended_x[appended_x %in% samples]

all_tested_genes <- scan('all_tested_genes.txt','character')

###################################################################################
# EXECUTE TEST #
#####################################################################################

perform_enrichment_test(data_table=transcriptomes,
                             columns=samples,
                             timepoints=transcriptome_descriptions[,2],
                             gaps_row=row_breaks,
                             gaps_col=col_breaks,
                             output_name=filename,
                             make_heatmaps=make_heatmaps)
