library(magrittr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicFeatures)
library(stringr)
library(dplyr)
library(data.table)
library("biomaRt")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(pals)
library(ggfortify)
library(ComplexHeatmap)

###Table from Count_Each_substitution

allAtoG<-read.table("completeOutput/a2g_sites/allAG_a2g.reditools_table.txt",sep="\t")
allAtoG[,7]<-as.numeric(as.character(str_extract(allAtoG[,7],"\\d+(?=,.\\d+])")))
allAtoG<-allAtoG[,c(1,2,4,7,5,15)]
colnames(allAtoG)<-c("chrom","position","strand","edits","coverage","Experiment")
allAtoG$strand <- case_when(
	allAtoG$strand == 1 ~ "+",
	allAtoG$strand == 0 ~ "-",
	TRUE ~ "*"
)

allAtoG$Edit_ratio <- allAtoG$edits / allAtoG$coverage
allAtoG$Location <- paste(allAtoG$chrom, allAtoG$position, allAtoG$strand)
write.table(allAtoG,sep="\t",quote=FALSE, row.names = FALSE,col.names = TRUE, "/completeOutput/a2g_sites/all_a2g_AG_readable.txt")

uniqueLocation<-table(allAtoG$Location)
atLeast4<-names(uniqueLocation[uniqueLocation>=3])
allAtoG<-allAtoG[allAtoG$Location %in% atLeast3,] 
write.table(allAtoG,sep="\t",quote=FALSE, row.names = FALSE,col.names = TRUE, "completeOutput/a2g_sites/all_a2g_AG_3sampEditreadable.txt")

#Table of all sites edit Ratio
allAtoG_editcast <- dcast(allAtoG, formula = Experiment ~  Location, value.var = "Edit_ratio", fill = 0, drop = FALSE)
row.names(allAtoG_editcast) <- allAtoG_editcast$Experiment
allAtoG_editcast$Experiment <- NULL
write.table(t(allAtoG_editcast),"completeOutput/a2g_sites/Human_a2g_editRatio.txt",sep = "\t",quote = FALSE,col.names=TRUE,row.names = TRUE)

#Table of all sites edit counts
allAtoG_editcast2 <- dcast(allAtoG, formula = Experiment ~  Location, value.var = "edits", fill = 0, drop = FALSE)
row.names(allAtoG_editcast2) <- allAtoG_editcast2$Experiment
allAtoG_editcast2$Experiment <- NULL
write.table(t(allAtoG_editcast2),"completeOutput/a2g_sites/Human_a2g_editCounts.txt",sep = "\t",quote = FALSE,col.names=TRUE,row.names = TRUE)

#### Bed of all edit Sites

#We need coverage for each site in each experiment, even the experiments that didn't have edits.

# The site locations can be split as if they were an input file.

allAtoG_editcast<-read.table("completeOutput/a2g_sites/Human_a2g_editRatio.txt",sep = "\t")
sites_split <- read.table(text = rownames(allAtoG_editcast), sep = " ", col.names = c("chrom", "position", "strand"), header = FALSE)

sites_split<-sites_split[grep("chr",sites_split$chrom),]
sites_granges <- GRanges(
	seqnames = sites_split$chrom,
	ranges = IRanges(
		start = as.numeric(sites_split$position),
		width = 1),
	strand = sites_split$strand,
	seqinfo = seqinfo(Hsapiens)
)
sites_granges <- sortSeqlevels(sites_granges)
sites_granges <- sort(sites_granges, ignore.strand = TRUE)
rtracklayer::export(sites_granges, "completeOutput/analysis/filtered_edit_sites.bed")

#### get depths for all samples at filtered edit sites
#### bash

posBAM=$(ls completeOutput/merged_alignmentsSPLIT/*/*positive.bam) 
negBAM=$(ls completeOutput/merged_alignmentsSPLIT/*/*negative.bam) 

module load samtools/1.9
printf "Chromosome\tPosition\tCoverage\tExperiment\n" > completeOutput/analysis/positiveSplit_edit_sites_coveragedepth.tsv
for bamfile in $(seq 0 59); do
	filename=$(basename ${posBAM[$bamfile]})
	experiment=${filename%%_*}
	samtools depth -b completeOutput/analysis/filtered_edit_sites.bed ${posBAM[$bamfile]} \
	| sed 's/$/\t'${experiment}'/g' \
	>> completeOutput/analysis/positiveSplit_edit_sites_coveragedepth.tsv
done

printf "Chromosome\tPosition\tCoverage\tExperiment\n" > completeOutput/analysis/negativeSplit_edit_sites_coveragedepth.tsv
for bamfile in $(seq 0 59); do
	filename=$(basename ${negBAM[$bamfile]})
	experiment=${filename%%_*}
	samtools depth -b completeOutput/analysis/filtered_edit_sites.bed ${negBAM[$bamfile]} \
	| sed 's/$/\t'${experiment}'/g' \
	>> completeOutput/analysis/negativeSplit_edit_sites_coveragedepth.tsv
done

#
posFulldf<-read.delim(paste0("completeOutput/analysis/negativeSplit_edit_sites_coveragedepth.tsv"))
	posFulldf$Location<-paste(posFulldf$Chromosome,posFulldf$Position)
	posFullcast<-dcast(posFulldf, formula = Experiment ~  Location, value.var = "Coverage", fill = 0, drop = FALSE)
	posFullcast<-data.frame(t(posFullcast))
	colnames(posFullcast)<-posFullcast[1,]
	posFullcast<-posFullcast[-1,]
	write.table(posFullcast,sep="\t",row.names=TRUE,col.names=TRUE,"completeOutput/analysis/filtered_edit_sites_negativeCount.tsv")
	
posFulldf<-read.delim(paste0("completeOutput/analysis/positiveSplit_edit_sites_coveragedepth.tsv"))
	posFulldf$Location<-paste(posFulldf$Chromosome,posFulldf$Position)
	posFullcast<-dcast(posFulldf, formula = Experiment ~  Location, value.var = "Coverage", fill = 0, drop = FALSE)
	posFullcast<-data.frame(t(posFullcast))
	colnames(posFullcast)<-posFullcast[1,]
	posFullcast<-posFullcast[-1,]
	write.table(posFullcast,sep="\t",row.names=TRUE,col.names=TRUE,"completeOutput/analysis/filtered_edit_sites_positiveCount.tsv")
	
	
#### Label GENE

sample<-NULL
allSites<-read.delim("completeOutput/a2g_sites/Human_a2g_editCount.txt")

temptbl<-data.frame(t(list2DF(strsplit(rownames(allSites)," ")))[,c(1,2,2,3)])
temptbl[,2]<-as.numeric(temptbl[,2])
temptbl[,3]<-as.numeric(temptbl[,3])
temptbl[,3]<-temptbl[,3]+1
temptbl[,1]<-gsub("chr","",temptbl[,1])
temptbl<-data.frame(temptbl,stringsAsFactors = FALSE)
temptbl<-transpose(as.list(temptbl))
all.genes<-NULL


for(i in 1:length(temptbl)){
  if(dim(genes[(genes$chromosome_name==temptbl[[i]][1]) & (genes$start_position<=as.numeric(temptbl[[i]][2]))&
               (genes$end_position>=as.numeric(temptbl[[i]][3])),])[1]==0){all.genes<-c(all.genes, list(c(paste0(temptbl[[i]][1]," ",temptbl[[i]][2]," ",temptbl[[i]][4]),"intergenic","intergenic")));next}
  

  primaryt<-genes[(genes$chromosome_name==temptbl[[i]][1]) & (genes$start_position<=as.numeric(temptbl[[i]][2]))&(genes$end_position>=as.numeric(temptbl[[i]][3])),c(3)][grep("principal",genes[(genes$chromosome_name==temptbl[[i]][1]) & (genes$start_position<=as.numeric(temptbl[[i]][2]))&(genes$end_position>=as.numeric(temptbl[[i]][3])),c(4)])[1]]
  all.genes<-c(all.genes,list(c(paste0(temptbl[[i]][1]," ",temptbl[[i]][2]," ",temptbl[[i]][4]),primaryt,genes[(genes$chromosome_name==temptbl[[i]][1]) & (genes$start_position<=as.numeric(temptbl[[i]][2]))&
                                                                                                              (genes$end_position>=as.numeric(temptbl[[i]][3])),c(2)])))
  if(i%%10000==0){print(paste0("At ",i))}
}
temp<-lapply(all.genes,function(x) x[1:3])
temp<-t(list2DF(temp))
temp[,1]<-paste0("chr",temp[,1])
write.table(temp,sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE,paste0("completeOutput/analysis/geneclassification.tsv"))


#### label repeat regions


sample<-NULL

genes<-read.delim("hg38_rmsk.bed") #from ucsc table browser
genes<-genes[,c(6:13)]
genes[,c(2:3)]<-genes[,c(2:3)]+1
temptbl<-data.frame(t(list2DF(strsplit(rownames(allAtoG_editcast)," ")))[,c(1,2,2,3)])
temptbl[,2]<-as.numeric(temptbl[,2])
temptbl[,3]<-as.numeric(temptbl[,3])
temptbl[,3]<-temptbl[,3]+1


temptbl<-data.frame(temptbl,stringsAsFactors = FALSE)
temptbl<-transpose(as.list(temptbl))
all.genes<-NULL
for(i in 1:length(temptbl)){
  if(dim(genes[(genes$genoName==temptbl[[i]][1]) & (genes$genoStart<=as.numeric(temptbl[[i]][2]))&
               (genes$genoEnd>=as.numeric(temptbl[[i]][2])),])[1]==0){all.genes<-c(all.genes, list(c(paste0(temptbl[[i]][1]," ",temptbl[[i]][2]," ",temptbl[[i]][4]),"None")));next}
  
  all.genes<-c(all.genes,list(c(paste0(temptbl[[i]][1]," ",temptbl[[i]][2]," ",temptbl[[i]][4]),
                                paste(genes[(genes$genoName==temptbl[[i]][1]) & (genes$genoStart<=as.numeric(temptbl[[i]][2]))&(genes$genoEnd>=as.numeric(temptbl[[i]][2])),c(7)],collapse="99"))))
   if(i%%5000==0){print(paste0("At ",i))}
}
write.table(t(list2DF(all.genes)),file = paste0("completeOutput/analysis/labeledEditSites_full.tsv"),sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#A table with each (primary) gene the edit falls into

#### label hyperedited regions

hypereditFiles<-dir("completeOutput/hyperedited_peaks/",".bed",full.names = TRUE)
hypers<-NULL
for(fl in 1:60){
fileHyper<-read.delim(hypereditFiles[fl],header = FALSE)
fileHyper[,2]<-fileHyper[,2]+1
fileHyper<-fileHyper[,c(1:3)]
fileHyper$sample<-str_extract(pattern = "(?<=//).+(?=.a2g)",hypereditFiles[fl])
hypers<-rbind(hypers,fileHyper)
}
hypers$site<-paste(hypers$V1,hypers$V2,hypers$V3)
uniquehyper<-hypers[!duplicated(hypers$site),]
write.table(uniquehyper,sep="\t",quote=FALSE, col.names=TRUE,row.names=FALSE,"completeOutput/analysis/unfileredHyper.tsv")

###xannotateHyper Rscript
#!/bin/bash
args = commandArgs(trailingOnly=TRUE)
chromo<-args[1]
allsites_info<-read.delim("completeOutput/analysis/primaryTranscript_locations_genomicSites.tsv")
allsites_info<-allsites_info[allsites_info$seqnames==chromo,]
hypers<-read.delim("completeOutput/analysis/unfileredHyper.tsv")
hypers<-hypers[hypers[,1]==chromo,]

hyperAnno<-apply(allsites_info[,],1,function(x) dim(table((hypers[,1] %in% x[[1]])&(x[[2]]>=hypers[,2])&(x[[2]]<=hypers[,3])))>1)
write.table(hyperAnno,paste0("/completeoutput/analysis/hyperAnno_",chromo,".tsv"),sep = "\t",col.names = TRUE, row.names = FALSE, quote = FALSE)


###


for chromo in chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr2 chr20 chr21 chr22 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrX chrY;
do sbatch xannotateHyper $chromo; done
allsites_info<-read.delim("completeOutput/analysis/primaryTranscript_locations_genomicSites.tsv")
allsites_info$hyperedited<-FALSE
for(chr in paste0("chr",c(1,10:19,2,20,21,22,3:9,"X","Y"))){
	hyperdf<-read.delim(paste0("completeOutput/analysis/hyperedits/hyperAnno_",chr,".tsv"),header=TRUE)
	allsites_info$hyperedited[allsites_info[,1]==chr]<-hyperdf[,1]
	}


write.table(rownames(allsites_info[allsites_info$hyperedited==TRUE,]),paste0("completeOutput/analysis/hyperedits/hypereditedSites.tsv"),sep = "\t",col.names = TRUE, row.names = FALSE, quote = FALSE)


#### label genetic regions

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl") #,host="apr2020.archive.ensembl.org"
genes = getBM(c("ensembl_gene_id","external_gene_name","ensembl_transcript_id","transcript_appris","chromosome_name","start_position","end_position","strand"), filters =c("chromosome_name","with_entrezgene"),values=list(c(1:22,"X","Y","MT"),TRUE), mart=mart)
genes<-genes[grep('principal',genes$transcript_appris),]
genes<-genes[!duplicated(genes$ensembl_gene_id),]

#primary transcripts

sites_split <- read.table(text = rownames(allAtoG), sep = " ", col.names = c("chrom", "position", "strand"), header = FALSE)

allsites_granges<-read.delim("completeOutput/analysis/filtered_edit_sites.bed",header=FALSE)

# Classify the sites
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
matching_exons <- cdsByOverlaps(txdb, allsites_granges, columns = c("GENEID","TXNAME"))
allsites_granges$is_exon <- overlapsAny(allsites_granges, matching_exons)
fiveUTR<-fiveUTRsByTranscript(txdb,use.names=TRUE)[gsub("\\.\\d+","",names(fiveUTRsByTranscript(txdb,use.names=TRUE))) %in% genes$ensembl_transcript_id]
allsites_granges$is_fiveUTR <- overlapsAny(allsites_granges, fiveUTR)
threeUTR<-threeUTRsByTranscript(txdb,use.names=TRUE)[gsub("\\.\\d+","",names(threeUTRsByTranscript(txdb,use.names=TRUE))) %in% genes$ensembl_transcript_id]
allsites_granges$is_threeUTR <- overlapsAny(allsites_granges, threeUTR)
intron<-intronsByTranscript(txdb,use.names=TRUE)[gsub("\\.\\d+","",names(intronsByTranscript(txdb,use.names=TRUE))) %in% genes$ensembl_transcript_id]
allsites_granges$is_intron <- overlapsAny(allsites_granges, intron)



allsites_info<-as.data.frame(allsites_granges)
allsites_info$location<-paste(allsites_info$seqnames,allsites_info$start,allsites_info$strand)

rownames(allsites_info)<-allsites_info$location

write.table(allsites_info,"completeOutput/analysis/primaryTranscript_locations_genomicSites.tsv",sep="\t",row.names = TRUE,quote=FALSE)



#### #### Full annotation Table

allsites_info<-read.delim("primaryTranscript_locations_genomicSites.tsv")
genelist<-read.delim("completeOutput/analysis/geneclassification.tsv",header = FALSE)
hypers<-read.delim("completeOutput/analysis/hyperedits/hypereditedSites.tsv",header=TRUE)
repeatlist<-read.delim("completeOutput/analysis/labeledEditSites_full.tsv",header=FALSE)
editRatio<-read.delim("completeOutput/analysis/Human_a2g_editRatio.txt")
editCount<-read.delim("completeOutput/analysis/Human_a2g_editCount.txt")
neg<-read.delim("completeOutput/analysis/filtered_edit_sites_negativeCount.tsv",header = TRUE)
pos<-read.delim("completeOutput/analysis/filtered_edit_sites_positiveCount.tsv",header = TRUE)
 
 fulldf<-neg+pos
 fulldf<-fulldf[str_extract(allsites_info$location,pattern=".+(?= (\\+|-|\\*))"),]

allsites_info$gene<-genelist[match(rownames(allsites_info),genelist[,1]),3]
rownames(genelist)<-genelist[,1]

allsites_info$genomicRegion<-"intergenic"
  for(row in 1:dim(allsites_info)[1]){ #rewrite using apply
    tmpdf<-data.frame(matrix(0,ncol = 2,nrow = 16),stringsAsFactors = FALSE)
  colnames(tmpdf)<-c("genomicLocation","which")
  tmpdf[,1]<-c("only CDS","only introns","only 5'UTR","only 3'UTR","CDS and 5' UTR","CDS and 3' UTR","CDS and introns",
               "CDS and both UTRs","CDS, introns, and 5' UTR","CDS, introns, and 3' UTR","NOT CDS, UTRs, or introns","3'UTR and intron","5'UTR and intron","both UTR","is all","intron and both UTRs")
  tmpdf[1,2]<-dim(subset(allsites_info[row,], (is_exon & !is_fiveUTR & !is_threeUTR & !is_intron)))[1]
  tmpdf[2,2]<-dim(subset(allsites_info[row,], (!is_exon & !is_fiveUTR & !is_threeUTR & is_intron)))[1]
  tmpdf[3,2]<-dim(subset(allsites_info[row,], (!is_exon & is_fiveUTR & !is_threeUTR & !is_intron)))[1]
  tmpdf[4,2]<-dim(subset(allsites_info[row,], (!is_exon & !is_fiveUTR & is_threeUTR & !is_intron)))[1]
  tmpdf[5,2]<-dim(subset(allsites_info[row,], (is_exon & is_fiveUTR & !is_threeUTR & !is_intron)))[1]          
  tmpdf[6,2]<-dim(subset(allsites_info[row,], (is_exon & !is_fiveUTR & is_threeUTR & !is_intron)))[1]             
  tmpdf[7,2]<-dim(subset(allsites_info[row,], (is_exon & !is_fiveUTR & !is_threeUTR & is_intron)))[1]
  tmpdf[8,2]<-dim(subset(allsites_info[row,], (is_exon & is_fiveUTR & is_threeUTR & !is_intron)))[1]
  tmpdf[9,2]<-dim(subset(allsites_info[row,], (is_exon & is_fiveUTR & !is_threeUTR & is_intron)))[1]
  tmpdf[10,2]<-dim(subset(allsites_info[row,], (is_exon & !is_fiveUTR & is_threeUTR & is_intron)))[1]
  tmpdf[11,2]<-dim(subset(allsites_info[row,], (!is_exon & !is_fiveUTR & !is_threeUTR & !is_intron)))[1]  
  tmpdf[12,2]<-dim(subset(allsites_info[row,], (!is_exon & !is_fiveUTR & is_threeUTR & is_intron)))[1] 
  tmpdf[13,2]<-dim(subset(allsites_info[row,], (!is_exon & is_fiveUTR & !is_threeUTR & is_intron)))[1]
  tmpdf[14,2]<-dim(subset(allsites_info[row,], (!is_exon & is_fiveUTR & is_threeUTR & !is_intron)))[1]
  tmpdf[15,2]<-dim(subset(allsites_info[row,], (is_exon & is_fiveUTR & is_threeUTR & is_intron)))[1]
  tmpdf[16,2]<-dim(subset(allsites_info[row,], (!is_exon & is_fiveUTR & is_threeUTR & is_intron)))[1]
  
  
    allsites_info$genomicRegion[row]<-tmpdf[which(tmpdf$which>0),1]
    if(row%%5000==0){print(row)}
  }
allsites_info$repeatRegion<-repeatlist[match(rownames(allsites_info),repeatlist[,1]),2]
allsites_info$hyperEdits<-"FALSE"
allsites_info$hyperEdits[(allsites_info$location %in% hypers[,1])]<-"TRUE"
allsites_info<-cbind(allsites_info,editCount[rownames(allsites_info),])


write.table(allsites_info,"completeOutput/analysis/full_annotation_table.tsv",sep="\t",row.names = TRUE,col.names = TRUE,quote=FALSE)  #### MASTER FILE


########### Files for normalization estimates

#For Repeats
repeattbl<-read.delim("references/hg38_rmsk.bed",header = TRUE)
repeattbl$repClass<-gsub("_","",repeattbl$repClass)
repeattbl$repClass<-gsub("Simplerepeat","SimpleRepeat",repeattbl$repClass)
repeattbl$repClass<-gsub("Lowcomplexity","LowComplexity",repeattbl$repClass)
repeattbl$repClass<-gsub("\\?","ZZZ",repeattbl$repClass)
replengthDF<-data.frame(matrix(0,nrow = 20,ncol = 2))
replengthDF[,1]<-unique(repeattbl$repClass)
for(rep in 1:20){
  
  matched<-grep(paste0(unique(repeattbl$repClass)[rep],"$"),repeattbl$repClass)
  replengthDF[rep,2]<-sum(repeattbl[matched,8]-repeattbl[matched,7])
  
}
replengthDF<-rbind(replengthDF,c("None",(3117275501-sum(replengthDF$X2)))) #3117275501 as total length
replengthDF$percent<-(3117275501/as.numeric(replengthDF[,2]))
replengthDF[,1]<-gsub("ZZZ","?",replengthDF[,1]) 

write.table(replengthDF,"completeOutput/analysis/repeatsNormTable.tsv",sep = "\t",col.names = TRUE, row.names = FALSE,quote=FALSE)

#######Genomic#####
genomicbed<-read.delim("completeOutput/analysis/UTRexport.txt",header = FALSE) #fasta of regions from uscs table browser
genomicbed<-genomicbed[(!duplicated(gsub("\\-\\d+","",genomicbed$Transcript.name))),]
genomicbed<-genomicbed[grep(pattern = ">",genomicbed$V1),1]
genomicbed<-t(list2DF(strsplit(genomicbed[grep(pattern = ">",genomicbed$V1),1],"\\|")))
genomicbed<-data.frame(genomicbed)
colnames(genomicbed)<-c("gene","5UTRSTART","5UTREND","3UTRSTART","3UTREND")
genomicbed[,1]<-gsub(">","",genomicbed[,1])
genomicbed<-genomicbed[genomicbed$gene %in% genes$ensembl_gene_id,]
UTR3len<-0
UTR5len<-0

for(row in 1:dim(genomicbed)[1]){
  if((length(grep("E",genomicbed[row,c(4)]))!=0)|suppressWarnings(is.na(as.numeric(genomicbed[row,4])))){}else{
  if(length(grep(";",genomicbed[row,c(4)]))==0){UTR3len<-UTR3len+(as.numeric(genomicbed[row,5])-as.numeric(genomicbed[row,4]))}else{
    tempdf<-rbind(genomicbed[row,4],genomicbed[row,5])
    tempdf<-data.frame(list2DF(strsplit(tempdf,";")))
    UTR3len<-UTR3len+sum(as.numeric(tempdf[,2])-as.numeric(tempdf[,1]))
  }}
  if((length(grep("E",genomicbed[row,c(2)]))!=0)|suppressWarnings(is.na(as.numeric(genomicbed[row,2])))){}else{
  if(length(grep(";",genomicbed[row,c(2)]))==0){UTR5len<-UTR5len+(as.numeric(genomicbed[row,3])-as.numeric(genomicbed[row,2]))}else{
    tempdf<-rbind(genomicbed[row,2],genomicbed[row,3])
    tempdf<-data.frame(list2DF(strsplit(tempdf,";")))
    UTR5len<-UTR5len+sum(as.numeric(tempdf[,2])-as.numeric(tempdf[,1]))
  }}
  
  if(row%%5000==0){print(row)}
}
fulllength<-3117275501
threeSum<-73888431
fiveSum<-6582662
cdsSum<-round(3117275501*0.011)

intronSum<-round(3117275501*0.24)
intergenicSum<-fulllength-sum(intronSum,fiveSum,threeSum,cdsSum)

fulllength/intergenicSum
fulllength/threeSum
fulllength/fiveSum
fulllength/intronSum
fulllength/cdsSum

genomicRegLenDF<-data.frame(matrix(0,ncol = 3,nrow = 5))
genomicRegLenDF[,1]<-names(sort(table(allsites_info$genomicRegion),decreasing = TRUE))[1:5][c(4,3,5,1,2)]
genomicRegLenDF[,2]<-c(threeSum,cdsSum,fiveSum,intronSum,intergenicSum)
genomicRegLenDF[,3]<-fulllength/genomicRegLenDF[,2]
write.table(genomicRegLenDF,"genomRegNormTable.tsv",sep = "\t",col.names = TRUE, row.names = FALSE,quote=FALSE)





################################################################################################
################################################################################################
# Generating figures

allSites<-read.delim("completeOutput/analysis/full_annotation_table.tsv",header = TRUE)

threshAtoGsites<-rownames(editRatio)[which(rowSums((editRatio>=0.1)&(editRatio<=0.9))>2)]
allSites$repeatRegion<-gsub("Simple_repeat","SimpleRepeat",allSites$repeatRegion)
allSites$repeatRegion<-gsub("Low_complexity","LowComplexity",allSites$repeatRegion)
allSites$repeatRegion<-gsub("99","_",allSites$repeatRegion)
 
repeatsNorm<-read.delim("completeOutput/analysis/repeatsNormTable.tsv",header=TRUE)
genomeNorm<-read.delim("completeOutput/analysis/genomRegNormTable.tsv",header = TRUE)
neg<-read.delim("completeOutput/analysis/filtered_edit_sites_negativeCount.tsv",header = TRUE)
pos<-read.delim("completeOutput/analysis/filtered_edit_sites_positiveCount.tsv",header = TRUE)

conditiongroups<-list(c(paste0("CKD0",seq(1,3,1))),c(paste0("DKD0",c(1:9)),
paste0("DKD",seq(10,30,1))),c(paste0("AKI0",seq(1,9,1)),paste0("AKI",seq(10,13,1))),
c(paste0("HRT0",seq(1,9,1)),paste0("HRT",seq(10,14,1)))) 
names(conditiongroups)<-c("CKD","DKD","AKI","HRT")
toremove<-c("DKD06","DKD02")	
allSites<-allSites[threshAtoGsites,!(colnames(allSites) %in% toremove)] ##These two samples have issues, 06 is dup of 05, 02 has depth mismatch 

for(rgrp in 1:length(toremove)){
	conditiongroups[[grep(toremove[rgrp],conditiongroups)]]<-conditiongroups[[grep(toremove[rgrp],conditiongroups)]][-grep(toremove[rgrp],conditiongroups[[grep(toremove[rgrp],conditiongroups)]])]}
 
fulldf<-neg+pos


#load in function 
source("BarplotsBySubset.R")

#User inputs to be generated
subsetLocation<-read.delim("two_threshold_sites.tsv",header = FALSE)[,1] #unique position ids
subsetName<-"above10below90"
subsetOUTDIR<-"completeOutput/analysis/AllEdit/with_two_threshold/"

my_function(subsetLocation = subsetLocation,subsetName = subsetName,subsetOUTDIR = subsetOUTDIR)


####### 
##Edit Percentage Fishers
##Find which edit sites are significant between conditions



  listdf<-list(NULL)
 
  #for the 4 conditiongroups 
  condList<-list(c(c(1,4)),c(c(2,4)),c(c(3,4)),c(c(1,2)),c(c(1,3)),c(c(2,3)),c(c(1,2,3,4))) #for "CKD","DKD","AKI","HRT"


narrowdf<-allSites
temp<-editRatio[allSites$location,colnames(allSites)[15:dim(allSites)[2]]]
temp<-mutate_all(temp, function(x) as.numeric(as.character(x)))
colnames(temp)<-colnames(allSites)[15:dim(allSites)[2]]
#All Sites where at least 3 samples in the comparison are edited over 10% but less than 90% and at least 3 samples have over 5 reads

  for(comp in 1:length(condList)){
	  if(comp==7){groupA<-c(conditiongroups[[condList[[comp]][1]]],conditiongroups[[condList[[comp]][2]]],conditiongroups[[condList[[comp]][3]]])
	  groupB<-c(conditiongroups[[condList[[comp]][4]]])}else{
	  groupA<-conditiongroups[[condList[[comp]][1]]]
	  groupB<-conditiongroups[[condList[[comp]][2]]]
		}
		print(comp)
    allSitesTemp<-cbind(narrowdf[,c(1:14)],narrowdf[,c(groupA,groupB)])
    allSitesTemp<-allSitesTemp[(rowSums(allSitesTemp[,15:dim(allSitesTemp)[2]]>5)>2),]
    editSites<-editRatio[rownames(allSitesTemp),c(groupA,groupB)]
    allSitesTemp<-allSitesTemp[(rowSums((editSites>0.1))>2)&(rowSums((editSites<0.9))>2),]

    genedf<-cbind(allSitesTemp[,c(1:14)],allSitesTemp[,c(groupA,groupB)])

    ddf<-data.frame(matrix(0,nrow = dim(genedf)[1],ncol = 5))
    if(comp==7){tempcondName<-"allvsHRT"}else{
    tempcondName<-paste0(names(conditiongroups)[condList[[comp]][1]],"_",names(conditiongroups)[condList[[comp]][2]])}
    colnames(ddf)<-c("location","gene",paste0("p.",tempcondName),paste0("p.adj.",tempcondName),paste0("LFC_",tempcondName))
   
    ddf[,1]<-genedf$location
    ddf[,2]<-genedf$gene
    
    editSites<-editRatio[genedf$location,colnames(allSitesTemp)[15:dim(allSitesTemp)[2]]]
      editSites<-cbind(round(rowMeans(editSites[,groupA])*100),round(rowMeans(editSites[,groupB])*100))
      editSites[editSites>100]<-100
      editSites[editSites<0]<-0
      
      
      counts_mat <- cbind(editSites[, 1], 100 - editSites[, 1], editSites[, 2], 100 - editSites[, 2])
      sample_index<-0
fisher_results <- apply(counts_mat, 1, function(row_counts) {
  p_value <- fisher.test(matrix(row_counts, nrow = 2))$p.value
  sample_index<-(sample_index+1)
   if (sample_index %% 5000 == 0) {
    message("At sample ", sample_index)
  }
  p_value
})

# Assign the results to the third column of the data frame
ddf[, 3] <- fisher_results
      
log_values <- data.frame(
  logA = ifelse(editSites[, 1] == 0, log2(0.0000001), log2(editSites[, 1])),
  logB = ifelse(editSites[, 2] == 0, log2(0.0000001), log2(editSites[, 2]))
)

# Calculate the difference between logA and logB
ddf[, 5] <- log_values$logA - log_values$logB

    ddf[,4]<-p.adjust(ddf[,3],method = "fdr")
    
    listdf[[comp]]<-ddf
    names(listdf)[comp]<-c("CKD_HRT","DKD_HRT","AKI_HRT","CKD_DKD","CKD_AKI","DKD_AKI","allvsHRT")[comp]
    allSitesTemp<-NULL
    some.genesTemp<-NULL
    totalCountTemp<-NULL
    }
  
######################################
##THIS PART USES EDIT PERCENT, GET EDIT PERCENT AS ITS OWN SUBSETDF (SUBSETeDF) AND MATCH UP COLUMNS or editPercent
#######################################

for(i in 1:length(condList)){
  listdf<-NULL
  for(reg in 1:length(unique(allsites_info$genomicRegion))){
    if(dim(editRatio[rownames(subsetDF[subsetDF$genomicRegion %in% unique(allsites_info$genomicRegion)[reg],]),])[1]==0){
      listdf[[reg]]<-0
    }else{
    temp<-editRatio[rownames(subsetDF[subsetDF$genomicRegion %in% unique(allsites_info$genomicRegion)[reg],]),][,condList[[i]]]
    temp<-mutate_all(temp, function(x) as.numeric(as.character(x)))
    temp[is.na(temp)]<-0
    listdf[[reg]]<-apply(temp,1,mean)}
  }
  names(listdf)<-unique(allsites_info$genomicRegion)
  longs<-NULL;for(b in 1:length(listdf)){longs<-append(longs,length(listdf[[b]]))}
  
  listdf<-listdf[longs>20]
  lengs<-NULL;for(b in 1:length(listdf)){lengs<-append(lengs,max(density(listdf[[b]])$y))};lengs<-max(lengs)
  
  pdf(height = 8,width = 8,useDingbats = FALSE, paste0(subsetOUTDIR,"/Histograms/By_GenomicRegion_PercentEditedSites_",condName[i],"_",subsetName,".pdf"))

  print(plot(density(listdf[[1]]),ylab = "Density",xlab="%Edit",xlim=c(0,1),main = paste0(condName[i]," % Edited in Site By Genomic Region ",subsetName),col=alpha(colorsused[1], 1),ylim=c(0,lengs*1.02)))
  print(lines(density(listdf[[2]]),col=alpha(colorsused[2], 1)))
  print(lines(density(listdf[[3]]),col=alpha(colorsused[3], 1)))
  if(length(listdf)>3){print(lines(density(listdf[[4]]),col=alpha(colorsused[4], 1)))}
  if(length(listdf)>4){print(lines(density(listdf[[5]]),col=alpha(colorsused[5], 1)))}
  if(length(listdf)>5){print(lines(density(listdf[[6]]),col=alpha(colorsused[6], 1)))}
  if(length(listdf)>6){ print(lines(density(listdf[[7]]),col=alpha(colorsused[7], 1)))}
  if(length(listdf)>7){ print(lines(density(listdf[[8]]),col=alpha(colorsused[8], 1)))}
  if(length(listdf)>8){  print(lines(density(listdf[[9]]),col=alpha(colorsused[9], 1)))}
  
  print(legend("topright",legend = names(listdf),fill=c(colorsused[1:length(listdf)])))
  dev.off()
}
dddf<-data.frame(matrix(0,nrow = dim(subsetDF)[1], ncol=11)) 

#### edit percentages by condition

condList<-conditiongroups
for(i in 1:4){
  
    condName<-names(conditiongroups)[i]
  band<-1
  tall<-10
  smalldf<-fulldf[str_extract(subsetDF$location, pattern=".+(?= (\\+|-|\\*))"),]


  dddf[,1]<-c(1:dim(subsetDF)[1])
  dddf[,(1+i)]<-rowSums(subsetDF[,condList[[i]]])/rowSums(smalldf[,condList[[i]]])
  matrixdf<-data.frame(matrix(0,nrow = dim(subsetDF)[1], ncol=length(conditiongroups[[i]]))) 
  matrixdf[,1]<-c(1:dim(subsetDF)[1])
  matrixdf[,2:dim(matrixdf)[2]]<-rowSums(subsetDF[,conditiongroups[[i]]])/rowSums(allsites_totalCounts[subsetDF$location,conditiongroups[[i]]])
  matrixdf[is.na(matrixdf)]<-0
  colnames(matrixdf)<-colnames(subsetDF[,conditiongroups[[i]]])
  pdf(height = 8,width = 8,useDingbats = FALSE, paste0(subsetOUTDIR,"/Histograms/Hist_PercentEditedSites_",condName[i],"_",subsetName,".pdf"))
  print(plot(density(matrixdf[,2]),ylab = "Density",xlab="%Edit",xlim=c(0,1),main = paste0(condName[i]," % Edited in Site ",subsetName),col=alpha(colorsused[1], 1),ylim=c(0,max(density(matrixdf[,2])[[2]],density(matrixdf[,3])[[2]],density(matrixdf[,4])[[2]],density(matrixdf[,5])[[2]],density(matrixdf[,6])[[2]])*1.1)))
  print(lines(density(matrixdf[,3]),col=alpha(colorsused[2], 1)))
  print(lines(density(matrixdf[,4]),col=alpha(colorsused[3], 1)))
  print(lines(density(matrixdf[,5]),col=alpha(colorsused[4], 1)))
  print(lines(density(matrixdf[,6]),col=alpha(colorsused[5], 1)))
  print(legend("topright",legend = c(conditiongroups[[1]]),fill=c(colorsused[1:6])))
  dev.off()
  
  pdf(height = 8,width = 8,useDingbats = FALSE, paste0(subsetOUTDIR,"/Histograms/Hist_PercentEditedSites_",condName[i],"_",subsetName,"_narrow.pdf"))
  print(plot(density(matrixdf[,2]),ylab = "Density",xlab="%Edit",xlim=c(0,band),main = paste0(condName[i]," % Edited in Site ",subsetName),col=alpha(colorsused[1], 1),ylim=c(0,tall)))
  print(lines(density(matrixdf[,3]),col=alpha(colorsused[2], 1)))
  print(lines(density(matrixdf[,4]),col=alpha(colorsused[3], 1)))
  print(legend("topright",legend = c(colnames(subsetDF)[c(9:38)][condList[[i]]]),fill=c(colorsused[1:3])))
  dev.off()
  write.table(matrixdf,sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE,paste0(subsetOUTDIR,"/Histograms/Hist_PercentEditedSites_",condName[i],"_",subsetName,".tsv"))
}

#### Azin1 figures

azinSites<-editRatio["chr8 102829408 -",]
 pdf(height = 8,width = 8,useDingbats = FALSE, paste0("completeOutput/analysis/PercentEditedSites_Azin1_chr8_102829408.pdf"))

  print(plot(density(rowMeans(azinSites[,grep("AKI",colnames(azinSites))]))
,ylab = "Density",xlab="%Edit",xlim=c(0,1),main = paste0("Human samples % Edited in Site By Condition"),col=alpha(colorsused[1], 1),ylim=c(0,13)))

  print(lines(density(rowMeans(azinSites[,grep("DKD",colnames(azinSites))])),col=alpha(colorsused[3], 1)))
  print(lines(density(rowMeans(azinSites[,grep("HRT",colnames(azinSites))])),col=alpha(colorsused[4], 1)))
  print(legend("topright",legend = c("AKI","CKD","DKD","HRT"),fill=c(colorsused[1:4])))
  dev.off()
  
  df<-data.frame(matrix(0,ncol=4,nrow=2))
  colnames(df)<-c("AKI","CKD","DKD","HRT"); rownames(df)<-c("chr8_102829408")
  df[1,1]<-rowMeans(azinSites[1,grep("AKI",colnames(azinSites))])
  df[1,2]<-rowMeans(azinSites[1,grep("CKD",colnames(azinSites))])
  df[1,3]<-rowMeans(azinSites[1,grep("DKD",colnames(azinSites))])
  df[1,4]<-rowMeans(azinSites[1,grep("HRT",colnames(azinSites))])
  
   pdf(height = 8,width = 8,useDingbats = FALSE, paste0("completeOutput/analysis/PercentEditedSites_Azin1_barplot.pdf"))
  barplot(as.matrix(df),beside=T,col=c(colorsused[1:2]))
  legend("topleft",c("chr8_102829408","chr8_102838824"),fill=c(colorsused[1:2]))
  dev.off()

ddf<-data.frame(matrix(0,ncol=4,nrow=60))
colnames(ddf)<-c("sample","condition","site","edit%")
ddf[,1]<-colnames(azinSites)
ddf[,2]<-c(rep("AKI",13),rep("CKD",3),rep("DKD",30),rep("HRT",14))
ddf[,3]<-rownames(azinSites)[1]
ddf[,4]<-as.numeric(azinSites[1,])
ddf$condition<-factor(ddf$condition,labels=c("AKI","CKD","DKD","HRT"))
ddf$sample<-factor(ddf$sample,ddf$sample)
box_plot<-ggplot(ddf,aes(x="condition",y="edit%"))+geom_boxplot()


ddf<-data.frame(matrix(0,ncol=3,nrow=60))
colnames(ddf)<-c("condition","editPercent","samples")
ddf[,1]<-c(rep("AKI",13),rep("CKD",3),rep("DKD",30),rep("HRT",14))
ddf[,2]<-as.numeric(azinSites[1,])
ddf[,3]<-colnames(azinSites)


pdf(height = 8,width = 8,useDingbats = FALSE, paste0("completeOutput/analysis/PercentEdited_chr8_102829408_Azin1_boxplot.pdf"))
ggplot(ddf[-22,], aes(x = condition, y = editPercent, fill = condition)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(color = "black", size = 0.7, alpha = 0.9) +
  geom_text(aes(label = samples), size = 3, position = position_nudge(y = -0.02))
  dev.off()

ddf<-data.frame(matrix(0,ncol=2,nrow=60))
colnames(ddf)<-c("condition","editPercent")
ddf[,1]<-c(rep("AKI",13),rep("CKD",3),rep("DKD",30),rep("HRT",14))
ddf[,2]<-as.numeric(azinSites[2,])


pdf(height = 8,width = 8,useDingbats = FALSE, paste0("completeOutput/analysis/PercentEdited_chr8_102838824_Azin1_boxplot.pdf"))
ggplot(ddf,aes(x=condition,y=editPercent,fill=condition))+geom_boxplot(outlier.color=NA)+
    geom_jitter(color="black", size=0.7, alpha=0.9)
dev.off()



###############
# Generate volcano plots of sigdiff edit sites by condition

# Generate Heatmaps of sigDiff edit sites by condition, sorted by GFR (for samples with GFR)
###############

# get Group Info
clinicalTable<-read.csv("metadata.csv")
rownames(clinicalTable)<-clinicalTable$study_id
sampleTable<-read.delim("samples.tsv")
sampleTable$known.condition<-gsub("H-CKD","CKD",sampleTable$known.condition)
sampleTable$CondA<-sampleTable$known.condition
sampleTable$CondB<-sampleTable$SEX
sampleTable$CondC<-paste0(sampleTable$known.condition,"_",sampleTable$SEX)
clinicalTable<-clinicalTable[sampleTable[grep("\\W",sampleTable$KPMP.Study.ID),2],]
clinicalTable<-clinicalTable[!is.na(clinicalTable$lab_egfr_most_recentC),c(1,2,25,5,6)]
tmpt<-sampleTable[grep("\\W",sampleTable$KPMP.Study.ID),c(2,6)]
tmpt<-tmpt[(tmpt[,1] %in% clinicalTable[,1]),]
clinicalTable$Merged.BAM<-tmpt$Merged.BAM
clinicalTable$GFR<-"normal"
clinicalTable$GFR[clinicalTable$lab_egfr_most_recentC<60]<-"CKD3a"
clinicalTable$GFR[clinicalTable$lab_egfr_most_recentC<=45]<-"CKD3b"
clinicalTable$GFR[clinicalTable$lab_egfr_most_recentC<=30]<-"CKD4"
clinicalTable$GFR[clinicalTable$lab_egfr_most_recentC<15]<-"CKD5"

sampleTable$egfr<-NULL; sampleTable$GFR<-NULL; sampleTable$np_age<-NULL; sampleTable$raceC<-NULL; 
for(merged in 1:dim(clinicalTable)[1]){sampleTable$egfr[sampleTable$Merged.BAM==clinicalTable$Merged.BAM[merged]]<-clinicalTable$lab_egfr_most_recentC[merged]; 
	sampleTable$GFR[sampleTable$Merged.BAM==clinicalTable$Merged.BAM[merged]]<-clinicalTable$GFR[merged];
	sampleTable$np_age[sampleTable$Merged.BAM==clinicalTable$Merged.BAM[merged]]<-clinicalTable$np_age[merged];
	sampleTable$raceC[sampleTable$Merged.BAM==clinicalTable$Merged.BAM[merged]]<-clinicalTable$raceC[merged]
	}

sampleTable$AZIN1_chr8_102829408_editCount<-NULL; sampleTable$AZIN1_chr8_102829408_editRatio<-NULL;	
for(merged in 1:dim(sampleTable)[1]){
	sampleTable$AZIN1_chr8_102829408_editCount[merged]<-editCount["chr8 102829408 -",sampleTable$Merged.BAM[merged]]
	sampleTable$AZIN1_chr8_102829408_editRatio[merged]<-round(editRatio["chr8 102829408 -",sampleTable$Merged.BAM[merged]],4)
}
	
write.table(sampleTable[,c(1,2,3,4,5,6,7,11,12,13,17:22)],sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE,"Clinical.txt")

#####


dirpath<-"completeOutput/analysis/figures/SigDiffCond_thresh_LFCpercent"
clinicalTable<-read.delim("Clinical.txt")
tempclin<-clinicalTable[!is.na(clinicalTable$egfr),]
tempclin<-tempclin[order(tempclin$egfr,decreasing=TRUE),]
tempclin$Merged.BAM<-factor(tempclin$Merged.BAM,levels=tempclin$Merged.BAM)

pdf(paste0("completeOutput/analysis/figures/SigDiffCond_thresh_LFCpercent/egfrplot.pdf"),width = 12,height = 8,useDingbats = FALSE) 
 
ggplot(data = tempclin, aes(x = Merged.BAM, y = egfr,group=1)) +
  geom_path() +
  labs(x = "Samples", y = "egfr") +
  ggtitle("egfr value of samples") +
  theme_minimal()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_rect(aes(xmin = as.numeric(Merged.BAM) - 0.5, xmax = as.numeric(Merged.BAM) + 0.5, ymin = -Inf, ymax = Inf, fill = GFR), alpha = 0.2)
dev.off()
azinSite<-c(as.character("Azin1_102829408"),as.numeric(editRatio["chr8 102829408 -",]))
azinSite<-rep(0,61)
for(condname in conds){
  subsetLocation<-read.delim(condsfls[grep(condname,condsfls)])
  table <- cbind(paste(subsetLocation[,1],subsetLocation[,3],subsetLocation[,4]),subsetLocation[,c(5,6)])
colnames(table)[1]<-"position"
 rownames(table)<-table$position
lfcthresh<-1

 pdf(paste0(dirpath,"/editSite_",condname,"_Volcanoplot_LFCthresh",lfcthresh,".pdf"),width = 12,height = 8,useDingbats = FALSE) 
  volcanoplot(table,lfcthresh=lfcthresh,sigthresh=0.05,textcx=0.4)
  dev.off()
  extremes<-c(table[order(table$LFC,decreasing=TRUE),1][1:250],table[order(table$LFC,decreasing=TRUE),1][(dim(table)[1]-249):dim(table)[1]])
  tmpRatio<-editRatio[extremes,]
  tmpLoc<-allSites[extremes,c(10,11,12,13)]
  tmpLoc$LFC<-table[rownames(tmpLoc),3]
  testdf<-data.frame(rbind(cbind(tmpLoc[order(tmpLoc$repeatRegion,tmpLoc$LFC),4],tmpRatio[tmpLoc[order(tmpLoc$repeatRegion,tmpLoc$LFC),1],]),azinSite))
  testdf<-testdf[c(dim(testdf)[1],1:(dim(testdf)[1]-1)),]
  rownames(testdf)[1]<-"Azin1_102829408"
colnames(testdf)[1]<-"repeatRegion"
testdf[1,1]<-"Azin1_102829408"
testdf[1,2:61]<-as.numeric(editRatio["chr8 102829408 -",])
 pdf(paste0(dirpath,"/Heatmap_Top500LFC_editSite_byRepeat_",condname,"_LFCthresh",lfcthresh,"_plusAzin1.pdf"),width = 12,height = 8,useDingbats = FALSE) 
print(Heatmap(as.matrix(testdf[,-1]),col=colorRampPalette(c("navy", "white", "firebrick3"))(n = 256),row_split =testdf[,1],cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize = 0.4),row_title_rot = 0,heatmap_legend_param = list(title = "Percent Editing"), column_title =paste0("Percent Editing of Top 500 LFC DE sites by repeat Region for ",condname)))
dev.off()
#heatmap with GFR order

gfrdf<-testdf[,c("repeatRegion",clinicalTable[order(clinicalTable$egfr,decreasing=TRUE),"Merged.BAM"])]
 pdf(paste0(dirpath,"/Heatmap_Top500LFC_editSite_byRepeat_",condname,"_LFCthresh",lfcthresh,"_by_GFR_HiLow.pdf"),width = 12,height = 8,useDingbats = FALSE) 
print(Heatmap(as.matrix(gfrdf[,-1]),col=colorRampPalette(c("navy", "white", "firebrick3"))(n = 256),row_split =gfrdf[,1],cluster_columns = FALSE, cluster_rows = FALSE,
row_names_gp = gpar(fontsize = 0.4),row_title_rot = 0,heatmap_legend_param = list(title = "Percent Editing"), column_title =paste0("Percent Editing of Top 500 LFC DE sites by repeat Region for ",condname,"by GFR")))
dev.off()
if(condname=="AKI_HRT"){
	gfrdf<-testdf[,c("repeatRegion",clinicalTable[order(clinicalTable$egfr,decreasing=TRUE),"Merged.BAM"][grep("AKI",clinicalTable[order(clinicalTable$egfr,decreasing=TRUE),"Merged.BAM"])])]
 pdf(paste0(dirpath,"/Heatmap_Top500LFC_editSite_byRepeat_",condname,"_LFCthresh",lfcthresh,"_AKI_by_GFR_HiLow.pdf"),width = 12,height = 8,useDingbats = FALSE) 
print(Heatmap(as.matrix(gfrdf[,-1]),col=colorRampPalette(c("navy", "white", "firebrick3"))(n = 256),row_split =gfrdf[,1],cluster_columns = FALSE, cluster_rows = FALSE,
row_names_gp = gpar(fontsize = 0.4),row_title_rot = 0,heatmap_legend_param = list(title = "Percent Editing"), column_title =paste0("Percent Editing of Top 500 LFC DE sites by repeat Region for ",condname,"by GFR_AKI only")))
dev.off()
	}
	if(condname=="DKD_HRT"){
	gfrdf<-testdf[,c("repeatRegion",clinicalTable[order(clinicalTable$egfr,decreasing=TRUE),"Merged.BAM"][grep("DKD",clinicalTable[order(clinicalTable$egfr,decreasing=TRUE),"Merged.BAM"])])]
 pdf(paste0(dirpath,"/Heatmap_Top500LFC_editSite_byRepeat_",condname,"_LFCthresh",lfcthresh,"_DKD_by_GFR_HiLow.pdf"),width = 12,height = 8,useDingbats = FALSE) 
print(Heatmap(as.matrix(gfrdf[,-1]),col=colorRampPalette(c("navy", "white", "firebrick3"))(n = 256),row_split =gfrdf[,1],cluster_columns = FALSE, cluster_rows = FALSE,
row_names_gp = gpar(fontsize = 0.4),row_title_rot = 0,heatmap_legend_param = list(title = "Percent Editing"), column_title =paste0("Percent Editing of Top 500 LFC DE sites by repeat Region for ",condname,"by GFR_DKD only")))
dev.off()
	}
}


########################################
##### finding junction regions between genes (MOUSE)

library("biomaRt")
library(stringr)
library(data.table)
mart = useMart("ensembl", dataset="mmusculus_gene_ensembl",host="apr2020.archive.ensembl.org")
genes = getBM(c("ensembl_gene_id","mgi_symbol","ensembl_transcript_id","transcript_appris","chromosome_name","start_position","end_position","strand"), mart=mart)
uniqGenes<-unique(genes$mgi_symbol)
startEnd<-function(gene){
  temp<-genes[genes$mgi_symbol %in% gene,]
  if(temp$strand[1]=="-1"){c(gene,paste0(temp$chromosome_name[1],"_",max(temp$end_position),"_",min(temp$start_position),"_-"))}
  if(temp$strand[1]=="1"){c(gene,paste0(temp$chromosome_name[1],"_",min(temp$start_position),"_",max(temp$end_position),"_+"))}
}

startEndNormal<-function(gene){
  temp<-genes[genes$mgi_symbol %in% gene,]
  if(temp$strand[1]=="-1"){c(gene,paste0(temp$chromosome_name[1],"_",min(temp$start_position),"_",max(temp$end_position),"_-"))}else{c(gene,paste0(temp$chromosome_name[1],"_",min(temp$start_position),"_",max(temp$end_position),"_+"))}
}


generegion<-rapply(as.list(uniqGenes), startEndNormal)
generegionDF<-cbind(generegion[seq(1,length(generegion),2)],generegion[seq(2,length(generegion),2)])
generegionDF<-data.frame(generegionDF)
colnames(generegionDF)<-c("Gene","Chrm_Start_End_Strand")
write.table(generegionDF,sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE, "GeneStartEnd_normal.tsv")
#Genes with starts and ends and directions

geneStartEnd<-read.delim("GeneStartEnd_normal.tsv")

dfbed<-data.frame(matrix(0,nrow = dim(geneStartEnd)[1],ncol = 6))
dfbed[,1]<-paste0("chr",as.array(data.frame(t(list2DF((strsplit(geneStartEnd[,2],split="_")))))[,1]))
dfbed[,2]<-(as.numeric(data.frame(t(list2DF(strsplit(geneStartEnd[,2],split="_"))))[,2])-1)
dfbed[,3]<-as.numeric(data.frame(t(list2DF(strsplit(geneStartEnd[,2],split="_"))))[,3])
dfbed[,4]<-geneStartEnd[,1]
dfbed[,5]<-"*"
dfbed[,6]<-(data.frame(t(list2DF(strsplit(geneStartEnd[,2],split="_"))))[,4])
dfbed<-dfbed[!(dfbed$X1 %in% unique(dfbed$X1)[23:30]),]
dfbed<-dfbed[dfbed$X3-dfbed$X2>=1000,] #removes 36434 leaving 18657
write.table(dfbed[dfbed$X6 == "-",],sep = "\t",col.names = FALSE,row.names = FALSE,"GeneStartEnd_normal_neg.bed")
write.table(dfbed[dfbed$X6 == "+",],sep = "\t",col.names = FALSE,row.names = FALSE,"GeneStartEnd_normal_pos.bed")
write.table(dfbed,sep = "\t",col.names = FALSE,row.names = FALSE,"GeneStartEnd_normal_full.bed")


## BASH
module load bedtools
sort -k1,1 -k2,2n GeneStartEnd_normal_neg.bed > GeneStartEnd_normal_neg.sorted.bed

bedtools closest -io -D ref -a GeneStartEnd_normal_neg.sorted.bed -b GeneStartEnd_normal_neg.sorted.bed > GeneStartEnd_normal_neg_closest.bed



sort -k1,1 -k2,2n /N/project/jeredSlate/2022/GY_AtoI/outputs/annotations/GeneStartEnd_normal_pos.bed > /N/project/jeredSlate/2022/GY_AtoI/outputs/annotations/GeneStartEnd_normal_pos.sorted.bed
bedtools closest -io -D ref -a /N/project/jeredSlate/2022/GY_AtoI/outputs/annotations/GeneStartEnd_normal_pos.sorted.bed -b /N/project/jeredSlate/2022/GY_AtoI/outputs/annotations/GeneStartEnd_normal_pos.sorted.bed > /N/project/jeredSlate/2022/GY_AtoI/outputs/annotations/GeneStartEnd_normal_pos_closest.bed


module load bedtools
sort -k1,1 -k2,2n GeneStartEnd_normal_full.bed > GeneStartEnd_normal_full.sorted.bed
bedtools closest -io -D ref -a GeneStartEnd_normal_full.sorted.bed -b GeneStartEnd_normal_full.sorted.bed > GeneStartEnd_normal_full_closest.bed

##In R
negGTF<-read.delim("GeneStartEnd_normal_neg_closest.bed",header=FALSE)
negGTF<-negGTF[,c(1,1,1,3,8,1,12,1,10)]
colnames(negGTF)<-c("Chr","xx","geneTranscript","Start","End","dot","Strand","dot2","GeneID")
negGTF[,c(6,8)]<-"."
negGTF[,2]<-"HAVANA"
negGTF[,3]<-"gene"
negGTF[,9]<-paste0("gene_name \"",negGTF[,9],"\"")
if(){} #nelect ends/starts that are 40000 max
over40000<-(negGTF[,5]-negGTF[,4])>40000
over10000<-(negGTF[,5]-negGTF[,4])>10000
negGTF10000<-negGTF
negGTF10000[over10000,5]<-(negGTF10000[over10000,4]+10000)
negGTF40000<-negGTF
negGTF40000[over40000,5]<-(negGTF40000[over40000,4]+40000)
write.table(negGTF40000,sep="\t",quote=FALSE,col.names = FALSE,row.names = FALSE,"GeneStartEnd_reverseGenes_300ntPlus_40000max.gtf")
write.table(negGTF10000,sep="\t",quote=FALSE,col.names = FALSE,row.names = FALSE,"GeneStartEnd_reverseGenes_300ntPlus_10000max.gtf")

#See if edits exist in these genes or not
#compare to counts of these

library("biomaRt")
library(stringr)
library(data.table)
mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")
genes = getBM(c("ensembl_gene_id","mgi_symbol","chromosome_name","start_position","end_position","strand"), filters =
c("chromosome_name","with_entrezgene"),
values=list(c(1:22,"X","Y","MT"),TRUE), mart=mart)

testbedfullALL<-read.delim("GeneStartEnd_normal_full_closest.bed",header = FALSE)
allSites<-read.delim("AtoG_Edits_FullAnnotations_withComparisons.tsv",header = TRUE)
allSites$RepeatRegion<-gsub("Simple_repeat","SimpleRepeat",allSites$RepeatRegion)
allSites$RepeatRegion<-gsub("Low_complexity","LowComplexity",allSites$RepeatRegion)
allSites<-allSites[rowSums(allSites[,c(9:38)]>20)>2,]
threshSites<-allSites[rowSums(allSites[,c(9:38)]>20)>2,]

testbedfull<-testbedfullALL[testbedfullALL[,6]=='+'&testbedfullALL[,12]=='-'&testbedfullALL[,13]>=300&testbedfullALL[,13]<=40000,]
siteList<-NULL
for(site in 1:dim(threshSites)[1]){
  siteList<-rbind(siteList,dim(testbedfull[(testbedfull[,1]==threshSites[site,2])&(testbedfull[,3]<threshSites[site,3])&(testbedfull[,8]>threshSites[site,3]),])[1])
   
}

tempSites<-threshSites[which(siteList>0),]

testbedfull$edited<-"FALSE"
for(site in 1:dim(testbedfull)[1]){
testbedfull$edited[site]<-dim(tempSites[(tempSites[,2]==testbedfull[site,1])&(tempSites[,3]>testbedfull[site,3])&(tempSites[,3])<testbedfull[site,8],])[1]>0
}
testbedfull$edits<-"NA"
testbedfull$editedStrand<-"NA"
for(site in 1:dim(testbedfull)[1]){
  dfdim<-dim(tempSites[(tempSites[,2]==testbedfull[site,1])&(tempSites[,3]>testbedfull[site,3])&(tempSites[,3])<testbedfull[site,8],])[1]
  if(dfdim==0){next}
  testbedfull$edits[site]<-paste(tempSites[(tempSites[,2]==testbedfull[site,1])&(tempSites[,3]>testbedfull[site,3])&(tempSites[,3])<testbedfull[site,8],1],collapse = ";")
  testbedfull$editedStrand[site]<-str_extract("\\-|\\+",string=testbedfull$edits[site])[1]
  #site<-which(testbedfull$edited==TRUE)[18]
 
}
testbedfull$junctionReverse<-paste0(testbedfull[,10],"_",testbedfull[,4],"_",testbedfull[,12])
testbedfull$junctionForward<-paste0(testbedfull[,4],"_",testbedfull[,10],"_",testbedfull[,6])
write.table(testbedfull,"~/2023/facingJunctions.tsv",sep = "\t",quote=FALSE,row.names = FALSE,col.names = TRUE)

###For facing regions
#Regions where genes are facing eachother (+ and -) with an edit site and are between 300 and 40000nt in length

fullGTF<-testbedfull
#fullGTF<-fullGTF[fullGTF[,4] %in%genes$mgi_symbol[genes$strand==1],]
fullGTF<-fullGTF[,c(1,1,1,3,8,1,12,1,10)]
fullGTF[,9]<-paste0(testbedfull[,10],"_",testbedfull[,4],"_",testbedfull[,12])
colnames(fullGTF)<-c("Chr","xx","geneTranscript","Start","End","dot","Strand","dot2","GeneID")
fullGTF[,c(6,8)]<-"."
fullGTF[,2]<-"HAVANA"
fullGTF[,3]<-"gene"
fullGTF[,9]<-paste0("gene_name \"",fullGTF[,9],"\"")
#c(1,1,1,3,8,1,6,1,4)

fullGTF2<-testbedfull
fullGTF2<-fullGTF2[,c(1,1,1,3,8,1,6,1,4)]
fullGTF2[,9]<-paste0(testbedfull[,4],"_",testbedfull[,10],"_",testbedfull[,6])
colnames(fullGTF2)<-c("Chr","xx","geneTranscript","Start","End","dot","Strand","dot2","GeneID")
fullGTF2[,c(6,8)]<-"."
fullGTF2[,2]<-"HAVANA"
fullGTF2[,3]<-"gene"
fullGTF2[,9]<-paste0("gene_name \"",fullGTF2[,9],"\"")

fullGTF<-rbind(fullGTF,fullGTF2)
write.table(fullGTF,sep="\t",quote=FALSE,col.names = FALSE,row.names = FALSE,"facingREGIONS_300nt_40kmax.gtf")


##########same direction regions
#Regions where genes are facing the same direction(+ + or - -) with an edit site and are between 300 and 40000nt in length
allSites<-read.delim("AtoG_Edits_FullAnnotations_withComparisons.tsv",header = TRUE)
allSites$RepeatRegion<-gsub("Simple_repeat","SimpleRepeat",allSites$RepeatRegion)
allSites$RepeatRegion<-gsub("Low_complexity","LowComplexity",allSites$RepeatRegion)
allSites<-allSites[rowSums(allSites[,c(9:38)]>20)>2,]
threshSites<-allSites[rowSums(allSites[,c(9:38)]>20)>2,]
testbedfull<-testbedfullALL[((testbedfullALL[,6]=='+'&testbedfullALL[,12]=='+')|(testbedfullALL[,6]=='-'&testbedfullALL[,12]=='-'))&testbedfullALL[,13]>=300&testbedfullALL[,13]<=40000,]

siteList<-NULL
for(site in 1:dim(threshSites)[1]){
  siteList<-rbind(siteList,dim(testbedfull[(testbedfull[,1]==threshSites[site,2])&(testbedfull[,3]<threshSites[site,3])&(testbedfull[,8]>threshSites[site,3]),])[1])
   
}

tempSites<-threshSites[which(siteList>0),]

testbedfull$edited<-"FALSE"
for(site in 1:dim(testbedfull)[1]){
testbedfull$edited[site]<-dim(tempSites[(tempSites[,2]==testbedfull[site,1])&(tempSites[,3]>testbedfull[site,3])&(tempSites[,3])<testbedfull[site,8],])[1]>0
}
testbedfull$edits<-"NA"
testbedfull$editedStrand<-"NA"
for(site in 1:dim(testbedfull)[1]){
  dfdim<-dim(tempSites[(tempSites[,2]==testbedfull[site,1])&(tempSites[,3]>testbedfull[site,3])&(tempSites[,3])<testbedfull[site,8],])[1]
  if(dfdim==0){next}
  testbedfull$edits[site]<-paste(tempSites[(tempSites[,2]==testbedfull[site,1])&(tempSites[,3]>testbedfull[site,3])&(tempSites[,3])<testbedfull[site,8],1],collapse = ";")
  testbedfull$editedStrand[site]<-str_extract("\\-|\\+",string=testbedfull$edits[site])[1]
  #site<-which(testbedfull$edited==TRUE)[18]
}
testbedfull$junctionReverse<-paste0(testbedfull[,10],"_",testbedfull[,4],"_",testbedfull[,12])
testbedfull$junctionForward<-paste0(testbedfull[,4],"_",testbedfull[,10],"_",testbedfull[,6])
write.table(testbedfull,"~/2023/sameDirectionJunctions.tsv",sep = "\t",quote=FALSE,row.names = FALSE,col.names = TRUE)

fullGTF<-testbedfull
fullGTF<-fullGTF[,c(1,1,1,3,8,1,12,1,10)]
fullGTF[,9]<-paste0(testbedfull[,10],"_",testbedfull[,4],"_",testbedfull[,12])
colnames(fullGTF)<-c("Chr","xx","geneTranscript","Start","End","dot","Strand","dot2","GeneID")
fullGTF[,c(6,8)]<-"."
fullGTF[,2]<-"HAVANA"
fullGTF[,3]<-"gene"
fullGTF[,9]<-paste0("gene_name \"",fullGTF[,9],"\"")

write.table(fullGTF,sep="\t",quote=FALSE,col.names = FALSE,row.names = FALSE,"sameDirectionREGIONS_300nt_40kmax.gtf")


#### for intron and junction edit distribution

txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene::TxDb.Mmusculus.UCSC.mm10.ensGene
intron<-intronsByTranscript(txdb,use.names=TRUE)[names(intronsByTranscript(txdb,use.names=TRUE)) %in% alltranscripts$Transcript.stable.ID]
allintrons<-unlist(intron)

allintrons_df <- data.frame(
  seqnames = seqnames(allintrons),
  start = start(allintrons),
  end = end(allintrons),
  name = names(allintrons),
  other = "*",
  strand = strand(allintrons)
)


write.table(allintrons_df, file = "/allintrons.bed", sep = "\t", col.names = FALSE, row.names=FALSE, quote = FALSE)

##in bash
module load bedops/2.4.37
closest-features --closest thresholded_editSites.bed allintrons.bed > closestIntron.bed

#if closest intron has multiple results this will split those:
input_bed<-"closestIntron.bed"
output_bed<-"closestIntron_split.bed"

lines <- readLines(input_bed)

# Initialize an empty vector to store the split lines
split_lines <- character(0)

# Iterate through each line
for (line in lines) {
  # Check if the line contains a semicolon
  if (grepl(";", line)) {
    # Split the line by semicolon
    parts <- unlist(strsplit(line, ";"))
    
    # Iterate through the split parts
    for (part in parts) {
      # Create a copy of the original line and replace the last column with the current part
      if(part!=parts[1]){
            modified_line <- paste(prePart, part, sep = "")  # Append the part
  }else{
	  prePart<-part
	  prePart<-gsub("\\|.*$","|",prePart)
	  modified_line<-part}
      split_lines <- c(split_lines, modified_line)  # Append the modified line
    }
  } else {
    # If no semicolon is found, keep the original line unchanged
    split_lines <- c(split_lines, line)
  }
}

# Write the split lines to the output BED file
writeLines(split_lines, con = output_bed)

grep 'intron' whichIntron_all.bed > editIntrons.bed

#in R######

bed_file<-read.table("thresholded_editSites.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
all_intron<-read.table("allintrons.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

#which introns have edits which belong to thresholded list
whichIntrons <- read.delim("closestIntron.bed",stringsAsFactors=FALSE,header=FALSE)
intronPre<-whichIntrons[,c(1,3,6,7,8)]
intronPre[,3]<-gsub("\\|\\w+","",intronPre[,3])
# Want to figure out direction to find "start" and "end" distance

intronPre$Start<-NA
intronPre$End<-NA
for(site in 1:dim(intronPre)[1]){
	if(intronPre[site,3]=="+"){intronPre$Start[site]<-intronPre[site,4];intronPre$End[site]<-intronPre[site,5]}else{
		intronPre$Start[site]<-intronPre[site,5];intronPre$End[site]<-intronPre[site,4]
	}
}

editDistance<-data.frame((intronPre[,6] - intronPre[,2]),(intronPre[,7] - intronPre[,2]))


###Facing###
allSites<-read.delim("AtoG_Edits_FullAnnotations.tsv",header = TRUE)
editSame<-read.delim("~/2023/sameDirectionJunctions.tsv",header = TRUE)
editFacing<-read.delim("~/2023/facingJunctions.tsv",header = TRUE)

#edit sites in facing regions
editSitesSame<-unlist(str_split(editSame[editSame$edited==TRUE,]$edits,pattern = ";"))
table(allSites[allSites$location %in% editSitesSame,]$RepeatRegion)
editSitesFacing<-unlist(str_split(editFacing[editFacing$edited==TRUE,]$edits,pattern = ";"))
table(allSites[allSites$location %in% editSitesFacing,]$RepeatRegion)


#####


editSiteDistributionIntron<-function(editSiteList,RegionList,pdfname){

for(edit in 1:length(editSiteList)){
	editSite<-strsplit(editSiteList[edit],split=" ")[[1]]
	editSite[2]<-as.numeric(editSite[2])
	editNumber<-as.numeric(editSite[2])
	region<-which((RegionList$chr == editSite[1])&(RegionList$start<=editNumber)&(RegionList$end>=editNumber))[1]
	if(is.na(region)){next}

	if(editSite[3]=="+"){
		start<-(RegionList[region,2]+6)
		end<-(RegionList[region,3]-40)
		bins <- list()
		bin_width <- ceiling((end - start + 1) / 20)
		for (i in 1:20) {
		bin_start <- start + (i - 1) * bin_width
		bin_end <- min(start + i * bin_width - 1, end)  # Ensure the last bin doesn't exceed the range
		if(i==20){bin_end<-end}
		bins[[i]] <- c(bin_start, bin_end)
		}
		if(any((RegionList[region,2]):(RegionList[region,2]+5) %in% editNumber)&(RegionList[region,5]<1)){RegionList[region,5]<-(RegionList[region,5]+1);
			endsbins<-rbind(endsbins,c(paste(editSite[1],editSite[2],editSite[3]),"GURAGN"))}
		for(binN in 1:20){
		if(any(bins[[binN]][1]:bins[[binN]][2] %in% editNumber)&(RegionList[region,(5+binN)]<1)){RegionList[region,(5+binN)]<-(RegionList[region,(5+binN)]+1)}
	
		
	}
		if(any((RegionList[region,3]-39):(RegionList[region,3]-10) %in% editNumber)&(RegionList[region,26]<1)){RegionList[region,26]<-(RegionList[region,26]+1);
			endsbins<-rbind(endsbins,c(paste(editSite[1],editSite[2],editSite[3]),"end30"))}
		if(any((RegionList[region,3]-9):(RegionList[region,3]) %in% editNumber)&(RegionList[region,27]<1)){RegionList[region,27]<-(RegionList[region,27]+1);
			endsbins<-rbind(endsbins,c(paste(editSite[1],editSite[2],editSite[3]),"end10"))}
}else{
		start<-(RegionList[region,3]-6)
		end<-(RegionList[region,2]+40)
		bins <- list()
		bin_width <- ceiling((end - start + 1) / 20)
		for (i in 1:20) {
		bin_start <- start + (i - 1) * bin_width
		bin_end <- max(start + i * bin_width - 1, end)  # Ensure the last bin doesn't exceed the range
		if(i==20){bin_end<-end}
		bins[[i]] <- c(bin_start, bin_end)
		}
		if(any((RegionList[region,3]-5):(RegionList[region,3]) %in% editNumber)&(RegionList[region,5]<1)){RegionList[region,5]<-(RegionList[region,5]+1);
			endsbins<-rbind(endsbins,c(paste(editSite[1],editSite[2],editSite[3]),"GURAGN"))}
	for(binN in 1:20){	
		if(any(bins[[binN]][1]:bins[[binN]][2] %in% editNumber)&(RegionList[region,(5+binN)]<1)){RegionList[region,(5+binN)]<-(RegionList[region,(5+binN)]+1)}	
		
}
		if(any((RegionList[region,2]+39):(RegionList[region,2]+10) %in% editNumber)&(RegionList[region,26]<1)){RegionList[region,26]<-(RegionList[region,26]+1);
			endsbins<-rbind(endsbins,c(paste(editSite[1],editSite[2],editSite[3]),"end30"))}
		if(any((RegionList[region,2]):(RegionList[region,2]+9) %in% editNumber)&(RegionList[region,27]<1)){RegionList[region,27]<-(RegionList[region,27]+1);
			endsbins<-rbind(endsbins,c(paste(editSite[1],editSite[2],editSite[3]),"end10"))}
	}
}

totalBins<-colSums(RegionList[,5:27])
totalBins

pdf(paste0("JunctionRegions/",pdfname,".pdf"),width=12,height=12,useDingbats=FALSE)
plot(totalBins,type='l',xaxt='n',xlab="region Location",ylab="# edits",main=pdfname)
abline(v = seq(1.5,22.5,1), col = adjustcolor("grey", alpha = 0.3))
axis(1, at=1:23, labels = FALSE)
text(seq(1, 23, by=1), par("usr")[3] - 2, labels = names(totalBins), srt = 90, pos = 1, xpd = TRUE)

dev.off()
pdf(paste0("JunctionRegions/",pdfname,"_bar.pdf"),width=12,height=12,useDingbats=FALSE)

barplot(totalBins, beside = TRUE, col = "blue", 
        xlab = "region Location", ylab = "# edits", 
        main = pdfname, names.arg = names(totalBins),     
        las = 2,
        cex.names = 1)

dev.off()
}

#######################################

####For junctions

######################################


editSiteDistributionJunction<-function(editSiteList,RegionList,pdfname){
	
for(edit in 1:length(editSiteList)){
	editSite<-strsplit(editSiteList[edit],split=" ")[[1]]
	editSite[2]<-as.numeric(editSite[2])
	editNumber<-as.numeric(editSite[2])
	region<-which((RegionList$chr == editSite[1])&(RegionList$start<=editNumber)&(RegionList$end>=editNumber))[1]

	if(editSite[3]=="+"){
		start<-(RegionList[region,2])
		end<-(RegionList[region,3])
		bins <- list()
		bin_width <- ceiling((end - start + 1) / 20)
		for (i in 1:20) {
		bin_start <- start + (i - 1) * bin_width
		bin_end <- min(start + i * bin_width - 1, end)  # Ensure the last bin doesn't exceed the range
		if(i==20){bin_end<-end}
		bins[[i]] <- c(bin_start, bin_end)
		}
		
		for(binN in 1:20){
		if(any(bins[[binN]][1]:bins[[binN]][2] %in% editNumber)){RegionList[region,(4+binN)]<-(RegionList[region,(4+binN)]+1)}
	
	}
		
}else{
		start<-(RegionList[region,3])
		end<-(RegionList[region,2])
		bins <- list()
		bin_width <- ceiling((end - start + 1) / 20)
		for (i in 1:20) {
		bin_start <- start + (i - 1) * bin_width
		bin_end <- max(start + i * bin_width - 1, end)  # Ensure the last bin doesn't exceed the range
		if(i==20){bin_end<-end}
		bins[[i]] <- c(bin_start, bin_end)
		}

		
	for(binN in 1:20){	
		if(any(bins[[binN]][1]:bins[[binN]][2] %in% editNumber)){RegionList[region,(4+binN)]<-(RegionList[region,(4+binN)]+1)}	
		
}
	}
}

totalBins<-colSums(RegionList[,5:24])
totalBins

pdf(paste0("JunctionRegions/",pdfname,".pdf"),width=12,height=12,useDingbats=FALSE)
plot(totalBins,type='l',xaxt='n',xlab="region Location",ylab="# edits",main=pdfname)
abline(v = seq(1.5,19.5,1), col = adjustcolor("grey", alpha = 0.3))
axis(1, at=1:20, labels = FALSE)
text(seq(1, 20, by=1), par("usr")[3] - 2, labels = names(totalBins), srt = 90, pos = 1, xpd = TRUE)

dev.off()
pdf(paste0("JunctionRegions/",pdfname,"_bar.pdf"),width=12,height=12,useDingbats=FALSE)

barplot(totalBins, beside = TRUE, col = "blue", 
        xlab = "region Location", ylab = "# edits", 
        main = pdfname, names.arg = names(totalBins),     
        las = 2,
        cex.names = 1)

dev.off()
}


#########################################
## FACING
########################################


editSiteDistributionJunctionFacing<-function(editSiteList,RegionList,pdfname){
	
for(edit in 1:length(editSiteList)){
	editSite<-strsplit(editSiteList[edit],split=" ")[[1]]
	editSite[2]<-as.numeric(editSite[2])
	editNumber<-as.numeric(editSite[2])
	region<-which((RegionList$chr == editSite[1])&(RegionList$start<=editNumber)&(RegionList$end>=editNumber))[1]

		start<-(RegionList[region,2])
		end<-(RegionList[region,3])
		bins <- list()
		bin_width <- ceiling((end - start + 1) / 20)
		for (i in 1:20) {
		bin_start <- start + (i - 1) * bin_width
		bin_end <- min(start + i * bin_width - 1, end)  # Ensure the last bin doesn't exceed the range
		if(i==20){bin_end<-end}
		bins[[i]] <- c(bin_start, bin_end)
		}
		
		for(binN in 1:20){
		if(any(bins[[binN]][1]:bins[[binN]][2] %in% editNumber)){RegionList[region,(4+binN)]<-(RegionList[region,(4+binN)]+1)}
	
	}
		

}

totalBins<-colSums(RegionList[,5:24])
totalBins

pdf(paste0("JunctionRegions/",pdfname,".pdf"),width=12,height=12,useDingbats=FALSE)
plot(totalBins,type='l',xaxt='n',xlab="region Location",ylab="# edits",main=pdfname)
abline(v = seq(1.5,19.5,1), col = adjustcolor("grey", alpha = 0.3))
axis(1, at=1:20, labels = FALSE)
text(seq(1, 20, by=1), par("usr")[3] - 2, labels = names(totalBins), srt = 90, pos = 1, xpd = TRUE)

dev.off()
pdf(paste0("JunctionRegions/",pdfname,"_bar.pdf"),width=12,height=12,useDingbats=FALSE)

barplot(totalBins, beside = TRUE, col = "blue", 
        xlab = "region Location", ylab = "# edits", 
        main = pdfname, names.arg = names(totalBins),     
        las = 2,
        cex.names = 1)

dev.off()
}



smalldf<-data.frame(whichIntrons[,1],whichIntrons[,2],whichIntrons[,3],whichIntrons[,6])
colnames(smalldf)<-c("chr","start","end","strand")
smalldf<-smalldf[!is.na(smalldf[,2]),]
smalldf<-smalldf[!duplicated(smalldf[,c(1:4)]),]


smalldf<-cbind(smalldf,data.frame(matrix(0,nrow=dim(smalldf)[1],ncol=23)))
colnames(smalldf)[5:27]<-c("GURAGN",paste0("bin",1:20),"end30","end10")
smalldf<-smalldf[!duplicated(smalldf),]
RegionList<-smalldf
editSiteList<-siteList


endsbins<-data.frame(matrix(0,ncol=2,nrow=0))
colnames(endsbins)<-c("editSite","region")
pdfname<-"AllIntronEditLocations_onePerbin"
editSiteDistributionIntron(editSiteList=siteList,RegionList=smalldf,pdfname=pdfname)
