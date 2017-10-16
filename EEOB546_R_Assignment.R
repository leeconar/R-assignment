library(plyr)
library(ggplot2)
library(reshape2)
######### Part 1:

# Load the file
SNP_data<-read.table("fang_et_al_genotypes.txt", fill=TRUE, header = FALSE, sep="\t", stringsAsFactors = FALSE)
SNP_position<-read.table("snp_position.txt",fill=TRUE, header = TRUE, sep="\t", stringsAsFactors = FALSE)

# Select the ID, chromosome and starting position for SNP_position
SNP_position<-data.frame(SNP_position$SNP_ID,SNP_position$Chromosome,SNP_position$Position)
colnames(SNP_position)<-c("ID", "Chromosome","Position")

##### Condition 1: ZMMIL, ZMMLR, ZMMMR

# Select group ZMMIL, ZMMLR, ZMMMR
SNP_sele1<-as.data.frame(t(subset(SNP_data, SNP_data$V3=="Group"|SNP_data$V3=="ZMMIL"|SNP_data$V3=="ZMMLR"|SNP_data$V3=="ZMMMR")))
SNP_sele1<-SNP_sele1[4:nrow(SNP_sele1),]
colnames(SNP_sele1)[1]<-c("ID")

# Merge by IDs
sele1_merged<-merge(SNP_position,SNP_sele1,by="ID")

####### Part I, Question 1: 
# Generate 10 files (1 for each chromosome) with SNP position sorted increasly, plus missing positions

# Change the missing position to ?
sele1_merged$Position<-sub("^$","?",sele1_merged$Position)

# Sort by position increasingly then select by each chromosome
sele1_merged<-arrange(sele1_merged,Position)

# Genereate the files
sele1_Chr1<-subset(sele1_merged, sele1_merged$Chromosome==1 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr2<-subset(sele1_merged, sele1_merged$Chromosome==2 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr3<-subset(sele1_merged, sele1_merged$Chromosome==3 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr4<-subset(sele1_merged, sele1_merged$Chromosome==4 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr5<-subset(sele1_merged, sele1_merged$Chromosome==5 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr6<-subset(sele1_merged, sele1_merged$Chromosome==6 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr7<-subset(sele1_merged, sele1_merged$Chromosome==7 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr8<-subset(sele1_merged, sele1_merged$Chromosome==8 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr9<-subset(sele1_merged, sele1_merged$Chromosome==9 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr10<-subset(sele1_merged, sele1_merged$Chromosome==10 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")

# Save files
write.table(sele1_Chr1,"Condition_1_Chr_1.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr2,"Condition_1_Chr_2.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr3,"Condition_1_Chr_3.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr4,"Condition_1_Chr_4.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr5,"Condition_1_Chr_5.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr6,"Condition_1_Chr_6.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr7,"Condition_1_Chr_7.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr8,"Condition_1_Chr_8.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr9,"Condition_1_Chr_9.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr10,"Condition_1_Chr_10.txt",sep="\t",row.names=FALSE)

# Remerge the file and change the missing position from ? to -
sele1_merged_replaced<-merge(SNP_position,SNP_sele1,by="ID")
sele1_merged_replaced$Position<-sub("^$","-",sele1_merged$Position)

# Sort by position decreasingly then select by each chromosome
sele1_merged_replaced<-arrange(sele1_merged_replaced,desc(Position))

# Generate the files
sele1_Chr1_replaced<-subset(sele1_merged, sele1_merged_replaced$Chromosome==1 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr2_replaced<-subset(sele1_merged, sele1_merged$Chromosome==2 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr3_replaced<-subset(sele1_merged, sele1_merged$Chromosome==3 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr4_replaced<-subset(sele1_merged, sele1_merged$Chromosome==4 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr5_replaced<-subset(sele1_merged, sele1_merged$Chromosome==5 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr6_replaced<-subset(sele1_merged, sele1_merged$Chromosome==6 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr7_replaced<-subset(sele1_merged, sele1_merged$Chromosome==7 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr8_replaced<-subset(sele1_merged, sele1_merged$Chromosome==8 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr9_replaced<-subset(sele1_merged, sele1_merged$Chromosome==9 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")
sele1_Chr10_replaced<-subset(sele1_merged, sele1_merged$Chromosome==10 & sele1_merged$Position!="multiple" &sele1_merged$Position!="unknown")

# Save the files
write.table(sele1_Chr1_replaced,"Condition_1_Chr_1_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr2_replaced,"Condition_1_Chr_2_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr3_replaced,"Condition_1_Chr_3_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr4_replaced,"Condition_1_Chr_4_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr5_replaced,"Condition_1_Chr_5_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr6_replaced,"Condition_1_Chr_6_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr7_replaced,"Condition_1_Chr_7_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr8_replaced,"Condition_1_Chr_8_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr9_replaced,"Condition_1_Chr_9_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele1_Chr10_replaced,"Condition_1_Chr_10_replaced.txt",sep="\t",row.names=FALSE)


####### Condition 2: ZMPBA, ZMPIL, ZMPJA

# Select group ZMPBA, ZMPIL, ZMPJA
SNP_sele2<-as.data.frame(t(subset(SNP_data, SNP_data$V3=="Group"|SNP_data$V3=="ZMPBA"|SNP_data$V3=="ZMPIL"|SNP_data$V3=="ZMPJA")))
SNP_sele2<-SNP_sele2[4:nrow(SNP_sele2),]
colnames(SNP_sele2)[1]<-c("ID")

# Merge by IDs
sele2_merged<-merge(SNP_position,SNP_sele2,by="ID")

####### Repeats


# Generate 10 files (1 for each chromosome) with SNP position sorted increasly, plus missing positions

# Change the missing position to ?
sele2_merged$Position<-sub("^$","?",sele2_merged$Position)

# Sort by position increasingly then select by each chromosome
sele2_merged<-arrange(sele2_merged,Position)

# Genereate the files
sele2_Chr1<-subset(sele2_merged, sele2_merged$Chromosome==1 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr2<-subset(sele2_merged, sele2_merged$Chromosome==2 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr3<-subset(sele2_merged, sele2_merged$Chromosome==3 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr4<-subset(sele2_merged, sele2_merged$Chromosome==4 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr5<-subset(sele2_merged, sele2_merged$Chromosome==5 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr6<-subset(sele2_merged, sele2_merged$Chromosome==6 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr7<-subset(sele2_merged, sele2_merged$Chromosome==7 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr8<-subset(sele2_merged, sele2_merged$Chromosome==8 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr9<-subset(sele2_merged, sele2_merged$Chromosome==9 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr10<-subset(sele2_merged, sele2_merged$Chromosome==10 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")

# Save files
write.table(sele2_Chr1,"Condition_2_Chr_1.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr2,"Condition_2_Chr_2.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr3,"Condition_2_Chr_3.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr4,"Condition_2_Chr_4.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr5,"Condition_2_Chr_5.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr6,"Condition_2_Chr_6.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr7,"Condition_2_Chr_7.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr8,"Condition_2_Chr_8.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr9,"Condition_2_Chr_9.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr10,"Condition_2_Chr_10.txt",sep="\t",row.names=FALSE)

# Remerge the file and change the missing position from ? to -
sele2_merged_replaced<-merge(SNP_position,SNP_sele2,by="ID")
sele2_merged_replaced$Position<-sub("^$","-",sele2_merged$Position)

# Sort by position decreasingly then select by each chromosome
sele2_merged_replaced<-arrange(sele2_merged_replaced,desc(Position))

# Generate the files
sele2_Chr1_replaced<-subset(sele2_merged, sele2_merged_replaced$Chromosome==1 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr2_replaced<-subset(sele2_merged, sele2_merged$Chromosome==2 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr3_replaced<-subset(sele2_merged, sele2_merged$Chromosome==3 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr4_replaced<-subset(sele2_merged, sele2_merged$Chromosome==4 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr5_replaced<-subset(sele2_merged, sele2_merged$Chromosome==5 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr6_replaced<-subset(sele2_merged, sele2_merged$Chromosome==6 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr7_replaced<-subset(sele2_merged, sele2_merged$Chromosome==7 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr8_replaced<-subset(sele2_merged, sele2_merged$Chromosome==8 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr9_replaced<-subset(sele2_merged, sele2_merged$Chromosome==9 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")
sele2_Chr10_replaced<-subset(sele2_merged, sele2_merged$Chromosome==10 & sele2_merged$Position!="multiple" &sele2_merged$Position!="unknown")

# Save the files
write.table(sele2_Chr1_replaced,"Condition_2_Chr_1_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr2_replaced,"Condition_2_Chr_2_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr3_replaced,"Condition_2_Chr_3_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr4_replaced,"Condition_2_Chr_4_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr5_replaced,"Condition_2_Chr_5_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr6_replaced,"Condition_2_Chr_6_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr7_replaced,"Condition_2_Chr_7_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr8_replaced,"Condition_2_Chr_8_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr9_replaced,"Condition_2_Chr_9_replaced.txt",sep="\t",row.names=FALSE)
write.table(sele2_Chr10_replaced,"Condition_2_Chr_10_replaced.txt",sep="\t",row.names=FALSE)






####### Part 2: 

## Plot 1: SNPs per chromosome

# Here, I don't really know which of the 2 conditions I should use
# So I reloaded the data and merged all the files, regardless of the Group
snps<-as.data.frame(t(SNP_data))
snps<-snps[4:nrow(snps),]
colnames(snps)[1]<-c("ID")
snps<-merge(SNP_position,snps,by="ID")

# Remove multiple and unknown
snps<-subset(snps, snps$Chromosome!="unknown" & snps$Chromosome!="multiple")

# Count how many snps per chromosome
snp_ct<-count(snps, c("Chromosome"))
colnames(snp_ct)[2]<-c("freq")

# Plotting
pl1<-ggplot(data=snp_ct,aes(x=snp_ct$Chromosome,y=snp_ct$freq))+geom_col()+
  theme(panel.background = element_rect(fill="#ffffff"),panel.border = element_rect(fill = NA,colour = "black"))+
  xlab("Chromosome")+ylab("Number of snps")+ggtitle("Snps per chromosome")

ggsave("plot1.pdf", pl1, width=10)

## Plot 2: 

# Reload the data, now this time, with headers
SNP_data2<-read.table("fang_et_al_genotypes.txt", fill=TRUE, header = TRUE, sep="\t", stringsAsFactors = FALSE)

# A function wrriten to determin the heterzygosity of each nucleotide pair
homo_hetero<-function(x){
  if (x== "A/A"| x== "T/T"|x == "C/C"| x== "G/G"){
    return ("homozygous")
  }
  else if (x=="?/?"){
    return ("NA")
  }
  else{
    return("heterozygous")
  }
}

# Reshape the file by melting
SNP_melted<-melt(SNP_data2,id=c("Sample_ID", "Group","JG_OTU"))
# Create the column with information of heterzygosity
SNP_melted$h_h<-as.character(lapply(SNP_melted$value, homo_hetero))

# Count the number of each homozygous, heterozygous and missing site for each group
homo_hete_ct<-count(SNP_melted,c("Group", "h_h"))
colnames(homo_hete_ct)[3]<-c("freq")

# Plotting, plus some fancy features
plt2<-ggplot(data=homo_hete_ct,aes(x="",y=freq,fill=h_h))+
  geom_bar(stat = "identity", position=position_fill())+
  facet_grid(~Group)+xlab("")+ylab("")+
  theme(axis.text.x = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(strip.text.x = element_text(size = 5, face="bold"))+
  theme(panel.background = element_rect(fill="#ffffff"),panel.border = element_rect(fill = NA,colour = "black"))+
  ggtitle("Heterozygosity per group")+scale_fill_discrete(name = "")

# Save the plot  
ggsave("plot2.pdf",plt2, width=20)    
