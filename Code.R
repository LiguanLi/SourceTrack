# SourceTrack

#R codes for data analyses and plotting
#including ARG abundance based sourcetracking, indicator, and correlation analysis.
#NOTE: the R script below has been simplified for illurstration

###################
#data preparation #
###################
install.packages("ggplot2")
install.packages("reshape")
install.packages("labdsv")
install.packages("RVAideMemoire")
install.packages("pheatmap")
library(ggplot2)
library(reshape)
library(labdsv)
library(RVAideMemoire)
library(pheatmap)


source('./src/SourceTracker.r') #load SourceTracker

database_str<-read.delim(file="database_str.txt",sep="\t",head=F) 
colnames(database_str)<-c("seqID","subtype","type") #table of database structure with column names of "seqID","subtype","type"

sample_info<-read.delim(file="sample_info.txt",sep='\t',header=T) #table of sample ecotype information with column names of "eco_type" and "sampleID"
sample_info$eco_type<-factor(sample_info$eco_type,levels = c("HF","AF","WA","NT"))
sample_info<-sample_info[order(sample_info$eco_type),]
sample_info$sampleID<-factor(sample_info$sampleID,levels=sample_info$sampleID)

ARG_abund<-read.delim(file="ARG_abund.txt",sep='\t',header=T) # table of ARG sequence abundance with column names of all seqIDs and row names of all sampleIDs
ARG_abundmatrix<-as.matrix(ARG_abund)
class(ARG_abundmatrix)<-"numeric"
ARG_abundmatrix<-ARG_abundmatrix[sample_info$sampleID,]

###############
#SourceTracker#
###############
#two functions for data processing of sourcetracking results
maxname<-function(data_frame){
  envs<-c("HF","AF","WA","NT","UN")
  max_name<-envs[which.max(data_frame[2:6])]
  return(max_name)}

maxratio<-function(data_frame){
  max_ratio<-data_frame[2:6][which.max(data_frame[2:6])]
  return(max_ratio)}

#leave-one-out strategy
pred_list<-list()
proptab_list<-list()
for (i in 1:nrow(ARG_abundmatrix))
{
  st<-sourcetracker(ARG_abundmatrix[-i,], sample_info$eco_type[-i])
  st_predic1<- predict(st,ARG_abundmatrix[i,], alpha1=0.001, alpha2=0.001)
  st_predic2 <- as.data.frame(st_predic1[[2]])
  st_predic2$sampleID<-rownames(ARG_abundmatrix)[i]
  st_predic2<-merge(st_predic2,sample_info,by="sampleID")
  pred_list[[i]]<-st_predic1
  proptab_list[[i]]<-st_predic2
}

proptab<-do.call("rbind",proptab_list)
proptab$predicted<-apply(proptab,1,maxname)
proptab$matched<-proptab$predicted == proptab$eco_type
proptab$predictedratio<-apply(proptab,1,maxratio)

#plotting
proptab_T<-subset(proptab,matched%in%"TRUE")
class(proptab_T$predictedratio)<-"double"
plotdata_stprob <- aggregate(proptab_T$predictedratio,by = list(eco_type=proptab_T$eco_type), FUN = function(x) c(mean = mean(x), sd = sd(x),n = length(x)))
plotdata_stprob <- do.call(data.frame, plotdata_stprob)
colnames(plotdata_stprob) <- c("eco_type", "mean", "sd", "n")
plotdata_stprob$eco_type<-ordered(plotdata_stprob$eco_type,levels=factor(c("HF","AF","WA","NT")))

png("prob_barplot.png", res = 300,width = 10, height = 10, units = 'in')
limits <- aes(ymax = plotdata_stprob$mean + plotdata_stprob$sd,ymin = plotdata_stprob$mean - plotdata_stprob$sd)
ggplot(data = plotdata_stprob, aes(x = eco_type, y = mean, fill = eco_type))+
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(limits, width = 0.3,color="gray50") +
  theme(axis.text.x=element_blank(), 
  		axis.ticks.x=element_blank(),
  		axis.title.x=element_blank(),
  		legend.position="none",
  		panel.background = element_rect(fill='gray96', colour='gray'),
  		text=element_text(size=25,color = "black"))+
  scale_y_continuous(breaks=seq(0,1,by=0.1))+
  scale_x_discrete(labels=c("HF","AF","WA","NT"))+
  scale_fill_manual(values = c("HF"="#FC7D00", "AF"="#FFC000","WA"="#5DA810","NT"="#15D9BA"))+
  theme_bw()+
  ylab("Probability")+
  xlab("")
graphics.off()


################
#Indicator ARGs#
################
clust<-sample_info$eco_type
seq_ind<-indval(ARG_abundmatrix,clust)
seq_pval_ind<-cbind(seq_ind$indval,seq_ind$pval)
seq_pval_ind$seqID<-rownames(seq_pval_ind)
seq_pval_ind<-merge(seq_pval_ind,database_str,by= "seqID",all.x=T)

sel_indseq<-list()
selpval_indseq<-subset(seq_pval_ind,seq_pval<=0.01)
eco_type<-c("HF","AF","WA","NT")
for (i in 1:length(eco_type))
{
	selpval_indseq<-selpval_indseq[order(selpval_indseq[,i],decreasing = T),]
	sel_indseq_tab<-selpval_indseq[1:30,]
	sel_indseq_tab$ind_cat<-eco_type[i]
	sel_indseq[[eco_type[i]]]<-sel_indseq_tab
}
seq_ind_summary<-Reduce(function(...) merge(..., all=T), sel_indseq)
seq_ind_summary$ind_cat<-factor(seq_ind_summary$ind_cat,levels=c("HF","AF","WA","NT"))
seq_ind_summary<-seq_ind_summary[order(seq_ind_summary$ind_cat),]

#plotting
plotdata_ind<-ARG_abund[,seq_ind_summary$seqID]
png("Indcator_heatmep.png", res = 300,width = 10, height = 10, units = 'in')
pheatmap(plotdata_ind, color = colorRampPalette(c("#AFAFAF","#FF9300"))(50),
			cluster_rows=AFLSE, cluster_cols=AFLSE, fontsize=5, fontsize_row=5, cellwidth=15, 
			cellheight=1, legend=TRUE,border_color = NA, show_rownames = F,show_colnames = T)
graphics.off()


###########################################
#Correlation by Linear Regression Analysis#
###########################################
cor_abundmatrix<-as.matrix(cbind(ARG_abund,rowSums(ARG_abund)))
cor_abundmatrix<-rbind(cor_abundmatrix,c(NA),c(NA))

for (i in 1:3502)
{
  lm_cor<-lm(as.numeric(unlist(cor_abundmatrix[1:656,3503]))~as.numeric(unlist(cor_abundmatrix[1:656,i])))
  cor_abundmatrix[657,i]<-summary(lm_cor)$r.squared
  if (sum(cor_abundmatrix[1:656,i]!=0)>=10 & cor_abundmatrix[678,i]>=0.5) cor_abundmatrix[658,i]<-"selcor"
}
R2_seq<-cor_abundmatrix[,!is.na(cor_abundmatrix[658,])]
R2_seq<-data.frame(t(R2_seq))
R2_seq$seqID<-rownames(R2_seq)

#plotting
R2_seq<-R2_seq[order(R2_seq[,657],decreasing = F),]
R2seq_abund<-ARG_abund[,R2_seq$seqID[1:5]]
R2seq_abund<-as.data.frame(R2seq_abund)
R2seq_abund$sampleID<-rownames(R2seq_abund)
R2seq_abund<-merge(R2seq_abund,sample_info,by="sampleID")
plotdata_R2seq<-melt(R2seq_abund)
plotdata_R2seq$variable2 = as.numeric(plotdata_R2seq$variable) + 10
plotdata_R2seq$eco_type<-factor(plotdata_R2seq$eco_type,levels=c("HF","AF","WA","NT"))
plotdata_R2seq<-plotdata_R2seq[order(plotdata_R2seq$eco_type),]
plotdata_R2seq$sampleID<-factor(plotdata_R2seq$sampleID,levels = sample_info$sampleID)
plotdata_R2seq<-plotdata_R2seq[order(plotdata_R2seq$sampleID),]

png("R2_corrmap.png", res = 300,width = 10, height = 10, units = 'in',bg = "transparent")
ggplot(plotdata_R2seq) +
  geom_tile(aes(x=sampleID, y=variable2, fill=factor(eco_type),alpha=value))+
  theme(panel.background=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  scale_fill_manual(values=c("#FC7D00","#FFC000","#5DA810","#15D9BA"),name="eco_type") +
  ylim(c(0,18))+
  coord_polar(theta="x")
graphics.off()  

