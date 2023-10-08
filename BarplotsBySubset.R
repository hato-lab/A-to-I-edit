###Want to input
##1) output location
##2) site subset input only
library(ggplot2)
library(stringr)
#subsetFile<-"/N/project/jeredSlate/2022/GY_AtoI/outputs/figures/GoodFigures/AllEdit/with_two_threshold/two_threshold_sites.tsv"  ##INPUT
#subsetOUTDIR<-"/N/project/jeredSlate/2022/GY_AtoI/outputs/figures/GoodFigures/AllEdit/with_two_threshold/"
#subsetDescriptor<-"above10below90"

#subsetLocation<-read.delim(subsetFile,header = FALSE)[,1]
#subsetName<-subsetDescriptor






my_function<-function(subsetLocation,subsetName,subsetOUTDIR){
  # allSites<-read.delim("/N/project/jeredSlate/2022/GY_AtoI/outputs/annotations/AtoG_Edits_Master.tsv",header = TRUE)
  # repeatsNorm<-read.delim("/N/project/jeredSlate/2022/GY_AtoI/outputs/annotations/repeatsNormTable.tsv",header=TRUE)
  # genomeNorm<-read.delim("/N/project/jeredSlate/2022/GY_AtoI/outputs/annotations/genomRegNormTable.tsv",header = TRUE)
  colorsused<-c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )[1:22]
  
subsetName<-subsetName
tmptest<-data.frame(matrix(NA,ncol = 12,nrow = 30))
colnames(tmptest)<-unique(allSites[,c(39)])
rownames(tmptest)<-colnames(allSites)[9:38]
subsetDF<-allSites[match(subsetLocation,allSites$location),]
for(i in 1:30){
  tmptest[i,]<-table(subsetDF[subsetDF[,(8+i)]>0,c(39)])[colnames(tmptest)]
}
tmptest[is.na(tmptest)]<-0
tmptest<-tmptest[,1:5]

#### Genomic Raw Figures
timepoint <- rep(rownames(tmptest),each=dim(tmptest)[2])
genomicLocation<- rep(colnames(tmptest),dim(tmptest)[1])
value<-NULL
for(i in 1:30){
  value<-append(value,as.numeric(tmptest[i,]))
}
data <- data.frame(timepoint,genomicLocation,value)
data$timepoint<-factor(data$timepoint,levels = rownames(tmptest))

nums<-NULL;for(i in seq(1,29,2)){nums<-append(nums,seq((i*5)-4,(i*5)))}
dir.create(paste0(subsetOUTDIR,"/BarPlots"))
dir.create(paste0(subsetOUTDIR,"/BarPlots/GenomicRegion"))
dir.create(paste0(subsetOUTDIR,"/BarPlots/RepeatRegion"))
dir.create(paste0(subsetOUTDIR,"/Histograms"))
pdf(paste0(subsetOUTDIR,"/BarPlots/GenomicRegion/threshold_EditSites_By_GenomicRegion_C_by_",subsetName,"_raw_percent.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="fill", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded edit sites by raw Genomic Region in Cytoplasm_",subsetName)))
dev.off()
write.table(data,paste0(subsetOUTDIR,"/BarPlots/GenomicRegion/threshold_EditSites_By_GenomicRegion_C_N_by_",subsetName,"_raw_values.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

pdf(paste0(subsetOUTDIR,"/BarPlots/GenomicRegion/threshold_EditSites_By_GenomicRegion_C_by_",subsetName,"_raw_values.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="stack", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded edit sites by raw Genomic Region in Cytoplasm_",subsetName)))
dev.off()


nums<-NULL;for(i in seq(2,30,2)){nums<-append(nums,seq((i*5)-4,(i*5)))}

pdf(paste0(subsetOUTDIR,"/BarPlots/GenomicRegion/threshold_EditSites_By_GenomicRegion_N_by_",subsetName,"_raw_percent.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="fill", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded edit sites by raw Genomic Region in Nuclear_",subsetName)))
dev.off()

pdf(paste0(subsetOUTDIR,"/BarPlots/GenomicRegion/threshold_EditSites_By_GenomicRegion_N_by_",subsetName,"_raw_values.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="stack", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded edit sites by raw Genomic Region in Nuclear_",subsetName)))
dev.off()
#### Genomic Normalized Figures

tmptest[,1]<-tmptest[,1]*genomeNorm[1,3]
tmptest[,2]<-tmptest[,2]*genomeNorm[2,3]
tmptest[,3]<-tmptest[,3]*genomeNorm[3,3]
tmptest[,4]<-tmptest[,4]*genomeNorm[4,3]
tmptest[,5]<-tmptest[,5]*genomeNorm[5,3]
 #normalized genomic values

timepoint <- rep(rownames(tmptest),each=5)
genomicLocation<- rep(colnames(tmptest),30)
value<-NULL
for(i in 1:30){
  value<-append(value,as.numeric(tmptest[i,]))
}
data <- data.frame(timepoint,genomicLocation,value)
data$timepoint<-factor(data$timepoint,levels = rownames(tmptest))

nums<-NULL;for(i in seq(1,29,2)){nums<-append(nums,seq((i*5)-4,(i*5)))}

pdf(paste0(subsetOUTDIR,"/BarPlots/GenomicRegion/threshold_EditSites_By_GenomicRegion_C_by_",subsetName,"_NORMALIZED_percent.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="fill", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded edit sites by Normalized Genomic Region in Cytoplasm_",subsetName)))
dev.off()

pdf(paste0(subsetOUTDIR,"/BarPlots/GenomicRegion/threshold_EditSites_By_GenomicRegion_C_by_",subsetName,"_NORMALIZED_values.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="stack", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded edit sites by Normalized Genomic Region in Cytoplasm_",subsetName)))
dev.off()


nums<-NULL;for(i in seq(2,30,2)){nums<-append(nums,seq((i*5)-4,(i*5)))}

pdf(paste0(subsetOUTDIR,"/BarPlots/GenomicRegion/threshold_EditSites_By_GenomicRegion_N_by_",subsetName,"_NORMALIZED_percent.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="fill", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded edit sites by Normalized Genomic Region in Nuclear_",subsetName)))
dev.off()

pdf(paste0(subsetOUTDIR,"/BarPlots/GenomicRegion/threshold_EditSites_By_GenomicRegion_N_by_",subsetName,"_NORMALIZED_values.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="stack", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded edit sites by Normalized Genomic Region in Nuclear_",subsetName)))
dev.off()
write.table(data,paste0(subsetOUTDIR,"/BarPlots/GenomicRegion/threshold_EditSites_By_GenomicRegion_C_N_by_",subsetName,"_NORMALIZED_values.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
#### Repeat Raw Figures


onlyFirstRepeat<-unlist(lapply(strsplit(allSites$RepeatRegion,"_"),'[[',1))
tmptest<-data.frame(matrix(NA,ncol = length(unique(onlyFirstRepeat)),nrow = 30))
colnames(tmptest)<-unique(onlyFirstRepeat)
rownames(tmptest)<-colnames(allSites)[9:38]
subsetDF$RepeatRegion<-unlist(lapply(strsplit(subsetDF$RepeatRegion,"_"),'[[',1))
for(i in 1:30){
  tmptest[i,]<-as.numeric(table(subsetDF[subsetDF[,(8+i)]>0,c(40)])[colnames(tmptest)])
}
tmptest[is.na(tmptest)]<-0

timepoint <- rep(rownames(tmptest),each=dim(tmptest)[2])
genomicLocation<- rep(colnames(tmptest),dim(tmptest)[1])
value<-NULL
for(i in 1:30){
  value<-append(value,as.numeric(tmptest[i,]))
}
data <- data.frame(timepoint,genomicLocation,value)
data$timepoint<-factor(data$timepoint,levels = rownames(tmptest))

nums<-NULL;for(i in seq(1,29,2)){nums<-append(nums,seq((i*dim(tmptest)[2])-(dim(tmptest)[2]-1),(i*dim(tmptest)[2])))}
pdf(paste0(subsetOUTDIR,"/BarPlots/RepeatRegion/threshold_EditSites_By_RepeatRegion_C_by_",subsetName,"_raw_percent.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="fill", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded raw edit sites by Repeat Region in Cytoplasm_",subsetName)))
dev.off()

pdf(paste0(subsetOUTDIR,"/BarPlots/RepeatRegion/threshold_EditSites_By_RepeatRegion_C_by_",subsetName,"_raw_value.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="stack", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded raw edit sites by Repeat Region in Cytoplasm_",subsetName)))
dev.off()

nums<-NULL;for(i in seq(2,30,2)){nums<-append(nums,seq((i*dim(tmptest)[2])-(dim(tmptest)[2]-1),(i*dim(tmptest)[2])))}

pdf(paste0(subsetOUTDIR,"/BarPlots/RepeatRegion/threshold_EditSites_By_RepeatRegion_N_by_",subsetName,"_raw_percent.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="fill", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded raw edit sites by Repeat Region in Nuclear_",subsetName)))
dev.off()

pdf(paste0(subsetOUTDIR,"/BarPlots/RepeatRegion/threshold_EditSites_By_RepeatRegion_N_by_",subsetName,"_raw_value.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="stack", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded raw edit sites by Repeat Region in Nuclear_",subsetName)))
dev.off()
write.table(data,paste0(subsetOUTDIR,"/BarPlots/RepeatRegion/threshold_EditSites_By_RepeatRegion_C_N_by_",subsetName,"_raw_values.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
### Repeats Normalized figures
for(rep in 1:length(unique(onlyFirstRepeat))){
  tmptest[,rep]<-tmptest[,rep]*repeatsNorm[repeatsNorm[,1]==colnames(tmptest)[rep],3] }
  
timepoint <- rep(rownames(tmptest),each=dim(tmptest)[2])
genomicLocation<- rep(colnames(tmptest),dim(tmptest)[1])
value<-NULL
for(i in 1:30){
  value<-append(value,as.numeric(tmptest[i,]))
}
data <- data.frame(timepoint,genomicLocation,value)
data$timepoint<-factor(data$timepoint,levels = rownames(tmptest))

nums<-NULL;for(i in seq(1,29,2)){nums<-append(nums,seq((i*dim(tmptest)[2])-(dim(tmptest)[2]-1),(i*dim(tmptest)[2])))}

pdf(paste0(subsetOUTDIR,"/BarPlots/RepeatRegion/threshold_EditSites_By_RepeatRegion_C_by_",subsetName,"_NORMALIZED_percent.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="fill", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded edit sites by Repeat Region in Cytoplasm_",subsetName)))
dev.off()

pdf(paste0(subsetOUTDIR,"/BarPlots/RepeatRegion/threshold_EditSites_By_RepeatRegion_C_by_",subsetName,"_NORMALIZED_value.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="stack", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded edit sites by Repeat Region in Cytoplasm_",subsetName)))
dev.off()

nums<-NULL;for(i in seq(2,30,2)){nums<-append(nums,seq((i*dim(tmptest)[2])-(dim(tmptest)[2]-1),(i*dim(tmptest)[2])))}

pdf(paste0(subsetOUTDIR,"/BarPlots/RepeatRegion/threshold_EditSites_By_RepeatRegion_N_by_",subsetName,"_NORMALIZED_percent.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="fill", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded edit sites by Repeat Region in Nuclear_",subsetName)))
dev.off()

pdf(paste0(subsetOUTDIR,"/BarPlots/RepeatRegion/threshold_EditSites_By_RepeatRegion_N_by_",subsetName,"_NORMALIZED_value.pdf"),height=8,width = 8,useDingbats = FALSE)
print(ggplot(data[nums,], aes(fill=genomicLocation, y=value, x=timepoint)) + 
        geom_bar(position="stack", stat="identity")+
        scale_fill_manual(values=colorsused)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("thresholded edit sites by Repeat Region in Nuclear_",subsetName)))
dev.off()
write.table(data,paste0(subsetOUTDIR,"/BarPlots/RepeatRegion/threshold_EditSites_By_RepeatRegion_C_N_by_",subsetName,"_NORMALIZED_values.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


###For histograms
condList<-list(
  c(1,3,5),c(7,9,11),c(13,15,17),c(19,21,23),c(25,27,29),c(2,4,6),c(8,10,12),c(14,16,18),c(20,22,24),c(26,28,30))
condName<-c("C_0hr","C_1hr","C_4hr","C_16hr","C_27hr","N_0hr","N_1hr","N_4hr","N_16hr","N_27hr")


#naming<-"above10below90"
for(i in 1:10){
  listdff<-NULL
  for(reg in 1:length(unique(allSites$GenomicLocation))){
    if(dim(do.call(rbind.data.frame,strsplit(subsetDF[subsetDF$GenomicLocation %in% unique(allSites$GenomicLocation)[reg],84],split = ":")))[1]==0){
      listdff[[reg]]<-0
    }else{
    temp<-do.call(rbind.data.frame,strsplit(subsetDF[subsetDF$GenomicLocation %in% unique(allSites$GenomicLocation)[reg],84],split = ":"))[,condList[[i]]]
    temp<-mutate_all(temp, function(x) as.numeric(as.character(x)))
    temp[is.na(temp)]<-0
    listdff[[reg]]<-apply(temp,1,mean)}
  }
  names(listdff)<-unique(allSites$GenomicLocation)
  longs<-NULL;for(b in 1:length(listdff)){longs<-append(longs,length(listdff[[b]]))}
  
  if(dim(list2DF(listdff))[1]<100){
	  listdff<-listdff[longs>1]
  lengs<-NULL;for(b in 1:length(listdff)){lengs<-append(lengs,max(listdff[[b]]))};lengs<-max(lengs)*2
	  }else{
	  listdff<-listdff[longs>20]
  lengs<-NULL;for(b in 1:length(listdff)){lengs<-append(lengs,max(listdff[[b]]))};lengs<-max(lengs)}
  
  
  pdf(height = 8,width = 8,useDingbats = FALSE, paste0(subsetOUTDIR,"/Histograms/By_GenomicRegion_PercentEditedSites_",condName[i],"_",subsetName,".pdf"))
  #print(hist(matrixdf[,2],type="l",ylab = "%Edit",xlab="Edit Site",ylim = c(0,1),main = paste0(condName[i]," % Edited in Site")))
  print(plot(density(listdff[[1]]),ylab = "Density",xlab="%Edit",xlim=c(0,1),main = paste0(condName[i]," % Edited in Site By Genomic Region ",subsetName),col=alpha(colorsused[1], 1),ylim=c(0,lengs)))
   if(length(listdff)>1){print(lines(density(listdff[[2]]),col=alpha(colorsused[2], 1)))}
  if(length(listdff)>2){print(lines(density(listdff[[3]]),col=alpha(colorsused[3], 1)))}
  if(length(listdff)>3){print(lines(density(listdff[[4]]),col=alpha(colorsused[4], 1)))}
  if(length(listdff)>4){print(lines(density(listdff[[5]]),col=alpha(colorsused[5], 1)))}
  if(length(listdff)>5){print(lines(density(listdff[[6]]),col=alpha(colorsused[6], 1)))}
  if(length(listdff)>6){ print(lines(density(listdff[[7]]),col=alpha(colorsused[7], 1)))}
  if(length(listdff)>7){ print(lines(density(listdff[[8]]),col=alpha(colorsused[8], 1)))}
  if(length(listdff)>8){  print(lines(density(listdff[[9]]),col=alpha(colorsused[9], 1)))}
  
  print(legend("topright",legend = names(listdf),fill=c(colorsused[1:length(listdf)])))
  dev.off()
}
dddf<-data.frame(matrix(0,nrow = dim(subsetDF)[1], ncol=11)) 
for(i in 1:10){
  
    
  band<-1
  tall<-10
  smalldf<-fulldf[match(subsetDF$location,fulldf$position),]
  
  dddf[,1]<-c(1:dim(subsetDF)[1])
  dddf[,(1+i)]<-rowSums(subsetDF[,c(9:38)[condList[[i]]]])/rowSums(smalldf[,c(10:39)[condList[[i]]]])
  matrixdf<-data.frame(matrix(0,nrow = dim(subsetDF)[1], ncol=1+length(condList[[i]]))) 
  matrixdf[,1]<-c(1:dim(subsetDF)[1])
  matrixdf[,2:(1+length(condList[[i]]))]<-subsetDF[,c(9:38)[condList[[i]]]]/smalldf[,c(10:39)[condList[[i]]]]
  matrixdf[is.na(matrixdf)]<-0
  colnames(matrixdf)<-colnames(subsetDF[,c(9:38)[condList[[i]]]])
  pdf(height = 8,width = 8,useDingbats = FALSE, paste0(subsetOUTDIR,"/Histograms/Hist_PercentEditedSites_",condName[i],"_",subsetName,".pdf"))
  #print(hist(matrixdf[,2],type="l",ylab = "%Edit",xlab="Edit Site",ylim = c(0,1),main = paste0(condName[i]," % Edited in Site")))
  print(plot(density(matrixdf[,2]),ylab = "Density",xlab="%Edit",xlim=c(0,1),main = paste0(condName[i]," % Edited in Site ",subsetName),col=alpha(colorsused[1], 1),ylim=c(0,max(density(matrixdf[,2])[[2]],density(matrixdf[,3])[[2]],density(matrixdf[,4])[[2]])*1.1)))
  print(lines(density(matrixdf[,3]),col=alpha(colorsused[2], 1)))
  print(lines(density(matrixdf[,4]),col=alpha(colorsused[3], 1)))
  print(legend("topright",legend = c(colnames(subsetDF)[c(9:38)][condList[[i]]]),fill=c(colorsused[1:3])))
  dev.off()
  
  pdf(height = 8,width = 8,useDingbats = FALSE, paste0(subsetOUTDIR,"/Histograms/Hist_PercentEditedSites_",condName[i],"_",subsetName,"_narrow.pdf"))
  #pdf(height = 8,width = 8,useDingbats = FALSE, paste0("/N/project/jeredSlate/2022/GY_AtoI/outputs/figures/GoodFigures/AllEdit/with_two_threshold/two_thresh_Hist_PercentEditedSites_",condName[i],"_",naming,"_narrow.pdf"))
  #print(hist(matrixdf[,2],type="l",ylab = "%Edit",xlab="Edit Site",ylim = c(0,1),main = paste0(condName[i]," % Edited in Site")))
  print(plot(density(matrixdf[,2]),ylab = "Density",xlab="%Edit",xlim=c(0,band),main = paste0(condName[i]," % Edited in Site ",subsetName),col=alpha(colorsused[1], 1),ylim=c(0,tall)))
  print(lines(density(matrixdf[,3]),col=alpha(colorsused[2], 1)))
  print(lines(density(matrixdf[,4]),col=alpha(colorsused[3], 1)))
  print(legend("topright",legend = c(colnames(subsetDF)[c(9:38)][condList[[i]]]),fill=c(colorsused[1:3])))
  dev.off()
  write.table(matrixdf,sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE,paste0(subsetOUTDIR,"/Histograms/Hist_PercentEditedSites_",condName[i],"_",subsetName,".tsv"))
}

pdf(height = 8,width = 8,useDingbats = FALSE, paste0(paste0(subsetOUTDIR,"/Histograms/Hist_PercentEditedSites_",subsetName,"_C_All_averages.pdf")))
#print(hist(matrixdf[,2],type="l",ylab = "%Edit",xlab="Edit Site",ylim = c(0,1),main = paste0(condName[i]," % Edited in Site")))
plot(density(dddf[,2]),ylab = "Density",xlab="%Edit",xlim=c(0,1),main = paste0("% Edited in Sites in Cytoplasmic Samples ",subsetName),col=alpha(colorsused[1], 1))
lines(density(dddf[,3]),col=alpha(colorsused[2], 1))
lines(density(dddf[,4]),col=alpha(colorsused[3], 1))
lines(density(dddf[,5]),col=alpha(colorsused[4], 1))
lines(density(dddf[,6]),col=alpha(colorsused[5], 1))
legend("topright",legend = c("0hr","1hr","4hr","16hr","27hr"),fill=c(colorsused[1:5]))
dev.off()

pdf(height = 8,width = 8,useDingbats = FALSE, paste0(paste0(subsetOUTDIR,"/Histograms/Hist_PercentEditedSites_",subsetName,"_N_All_averages.pdf")))
plot(density(dddf[,7]),ylab = "Density",xlab="%Edit",xlim=c(0,1),main = paste0("% Edited in Sites in Nuclear Samples ",subsetName),col=alpha(colorsused[1], 1))
lines(density(dddf[,8]),col=alpha(colorsused[2], 1))
lines(density(dddf[,9]),col=alpha(colorsused[3], 1))
lines(density(dddf[,10]),col=alpha(colorsused[4], 1))
lines(density(dddf[,11]),col=alpha(colorsused[5], 1))
legend("topright",legend = c("0hr","1hr","4hr","16hr","27hr"),fill=c(colorsused[1:5]))
dev.off()
write.table(dddf,sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE,paste0(subsetOUTDIR,"/Histograms/Hist_PercentEditedSites_FullAverages_",subsetName,".tsv"))

}
