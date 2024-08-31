library(tidyverse)
library(wesanderson)
library(ggpubr)
#violin plot to display data quality
#metabolomics
metabolomics_qc<-read.csv("metabolomics_violin_plot.csv",header =TRUE)
p<-ggplot(metabolomics_qc,aes(x=mode,y=cv,fill=mode))+
  geom_violin(adjust=0.7)+
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize=0.7,binwidth = 0.4,show.legend=FALSE,fill="black")+
  scale_fill_manual(values=wes_palette("Zissou1"))+
  labs(title="Metabolomics", x=NULL,y = "Coefficient of variation(%)")+
  theme(legend.position="none")+
  scale_x_discrete(limits=c("pos", "neg"))
#lipidomics
lipidomics_qc<-read.csv("lipidomics_violin_plot.csv",header =TRUE)
q<-ggplot(lipidomics_qc,aes(x=mode,y=cv,fill=mode))+
  geom_violin(adjust=0.7)+
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize=0.7,binwidth = 0.4,show.legend=FALSE,fill="black")+
  scale_fill_manual(values=wes_palette("Zissou1"))+
  labs(title="Lipidomics", x=NULL,y = "Coefficient of variation(%)")+
  theme(legend.position="none")+
  scale_x_discrete(limits=c("pos", "neg"))
ggarrange(p, q, 
          labels = c("A", "B"),
          ncol = 2, nrow =1)
ggsave("QCdistribution.png", width = 8, height =4, units = c("in"), dpi = 300)
#Display compounds class
library(dplyr)
library(RColorBrewer)
compounds<-read.csv("analysis_compounds.csv",header=TRUE)
compounds_sum<-compounds %>%
  group_by(platform,class) %>%
  summarize(count=n())
write.csv(compounds_sum,row.names=FALSE,file="compounds_sum1.csv")
#Here, some operations are done in EXCEL to simplify the table.
compounds_sum<-read.csv("compounds_sum1.csv",header=TRUE)
x_p<-c("Carnitines","Amino acids","Bile acids","Organic acids","Others","LPC","LPE","PC","PE","SM","FFA","Cer","CL","Co","DG","PI","TG")
y_p<-c(22,16,8,7,20,40,18,49,43,23,48,22,8,3,16,12,124)
ggplot(compounds_sum,aes(x=class,y=count,fill=platform))+
  geom_bar(stat="identity")+
  scale_fill_manual(values =c("#1f78b4","#b2df8a"))+
  scale_x_discrete(limits=c("Carnitines","Amino acids","Bile acids","Organic acids","Others","LPC","LPE","PC","PE","SM","FFA","Cer","CL","Co","DG","PI","TG"))+
  theme(legend.position=c(0.5,0.9),legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x=NULL,y = "Count")+
  annotate(geom="text", x=x_p, y=y_p+2, label=y_p)
ggsave("CompoundClass1.png", width = 13, height = 6, units = c("in"), dpi = 300)

#PCA PCAdata.csv has substituted NA as 0.
library(factoextra)
library(wesanderson)
PCAdata<-read.csv("PCAdata1.csv",header=TRUE,row.names = 1)
PCAdata.active<-PCAdata[,2:480]
res.pca <- prcomp(PCAdata.active, scale = TRUE,)
fviz_eig(res.pca)
PCAdata$group<-factor(PCAdata$group,levels=c('Control','PFO5DoDA-low','PFO5DoDA-high'))
group <- as.factor(PCAdata$group)
fviz_pca_ind(res.pca,
             col.ind =group, # color by group
             geom ="point",
             palette=c("#66c2a5", "#fc8d62", "#8da0cb"), addEllipses=TRUE,ellipse.level=0.6,
             legend.title="Group",invisible="quali"
)
ggsave("PCAscore1.png", width = 7, height = 6, units = c("in"), dpi = 300)

#Differential compounds analysis
#Shapiro test is used to test normallity of the samples.
PCAdata<-read.csv("PCAdata1.csv",header=TRUE,row.names = 1)
hist(subset(PCAdata, group == "PFO5DoDA-high")$LPC.22_6.sn.1,
     main = "LPC 22:6 sn 1 distribution",
     xlab = "LPC 22:6 sn 1"
)
shapiro.test(subset(PCAdata, group == "PFO5DoDA-high")$FA.14.0.)
#Compound boxplot
ggplot(PCAdata) +
  aes(x = group, y = Butyrylcarnitine) +
  geom_boxplot() +
  theme_minimal()+
  scale_x_discrete(limits=c("Control","PFO5DoDA-low","PFO5DoDA-high"))
compound<-colnames(PCAdata)
compound <- compound[! compound %in% c('group')]
low_PCAdata <- PCAdata %>% 
  filter(group == "Control"|group == "PFO5DoDA-low")
high_PCAdata <- PCAdata %>% 
  filter(group == "Control"|group == "PFO5DoDA-high")
library(plyr)
controlgauss.p<-ldply(
  compound,
  function(colname){
    p_val=shapiro.test(subset(PCAdata, group == "Control")[[colname]])$p.value
    return(data.frame(colname=colname, p_value=p_val))
  }
)
controlgauss.p<-controlgauss.p %>%
  mutate(is.less.than0.05=p_value<0.05)
#Mann-Whitney U test PCAdata above is more appropriate for analysis.

low_pvalue <- ldply(
  compound,
  function(colname) {
    p_val = wilcox.test(low_PCAdata[[colname]] ~ low_PCAdata$group)$p.value
    return(data.frame(colname=colname, p_value=p_val))
  })

high_pvalue<- ldply(
  compound,
  function(colname) {
    p_val = wilcox.test(high_PCAdata[[colname]] ~ high_PCAdata$group)$p.value
    return(data.frame(colname=colname, p_value=p_val))
  })
#Compute FDR from p_values
low_FDR<-p.adjust(low_pvalue$p_value, method ="BH")
high_FDR<-p.adjust(high_pvalue$p_value, method ="BH")

#Compute fold change
compound_data<-read.csv("CompoundsQuantity1.csv",header = TRUE)
compound_class<-read.csv("analysis_compounds.csv",header=TRUE)
compound_data <- compound_data %>% 
  mutate(
    ctrl_mean = (Ctrl.1.1+Ctrl.1.10+Ctrl.1.11+Ctrl.1.12+Ctrl.1.3
                 +Ctrl.1.4+Ctrl.1.5+Ctrl.1.6+Ctrl.1.7+Ctrl.1.8+Ctrl.1.9)/11,
    low_mean = (PFO5DoDA.2.1+PFO5DoDA.2.11+PFO5DoDA.2.2+PFO5DoDA.2.3
                +PFO5DoDA.2.4+PFO5DoDA.2.5+PFO5DoDA.2.6+PFO5DoDA.2.7
                +PFO5DoDA.2.8+PFO5DoDA.2.9)/10,
    high_mean=(PFO5DoDA.3.1+PFO5DoDA.3.10+PFO5DoDA.3.11+PFO5DoDA.3.12
               +PFO5DoDA.3.2+PFO5DoDA.3.4+PFO5DoDA.3.5+PFO5DoDA.3.6
               +PFO5DoDA.3.7+PFO5DoDA.3.9)/10,
    low_fc=low_mean/ctrl_mean,
    high_fc=high_mean/ctrl_mean
  ) %>%
  inner_join(compound_class,by="compounds")
compound_data<-data.frame(compound_data,low_pvalue=low_pvalue$p_value,
                             high_pvalue=high_pvalue$p_value,
                             low_FDR=low_FDR,high_FDR=high_FDR)
write.csv(compound_data,row.names=FALSE,file="Analyze_result1.csv")

#volcano plot
library(EnhancedVolcano)
compound_data<-read.csv("Analyze_result1.csv",header=TRUE)
compound_data <- compound_data %>% 
  mutate(
    log2low_fc=log2(low_fc),
    log2high_fc=log2(high_fc)
  )
a<- EnhancedVolcano(compound_data,
                     lab =compound_data$Standard.name, 
                     x = 'log2low_fc',
                     y = 'low_pvalue',
                    axisLabSize = 10,
                    title = 'PFO5DoDA(low dose) vs Control',
                    titleLabSize = 12,
                    pCutoff = 0.05,
                    FCcutoff = 0.5,
                    pointSize = 0.6,
                    labSize = 2.5,
                    col=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3'),
                    colAlpha = 0.9,
                    legendPosition = 'none',
                    subtitle =NULL,
                    caption = bquote(~Log[2]~ "fold change cutoff, 0.5; p-value cutoff, 0.05"),
                    captionLabSize = 9
                    )
b<-EnhancedVolcano(compound_data,
                   lab =compound_data$Standard.name, 
                   x = 'log2high_fc',
                   y = 'high_pvalue',
                   axisLabSize = 10,
                   title = 'PFO5DoDA(high dose) vs Control',
                   titleLabSize = 12,
                   pCutoff = 0.05,
                   FCcutoff = 0.5,
                   pointSize = 0.6,
                   labSize = 2.5,
                   col=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3'),
                   colAlpha = 0.9,
                   legendPosition = 'none',
                   subtitle =NULL,
                   caption = bquote(~Log[2]~ "fold change cutoff, 0.5; p-value cutoff, 0.05"),
                   captionLabSize = 9
)
ggarrange(a, b, 
          labels = c("A", "B"),
          ncol = 2, nrow =1)
ggsave("Volcano plot1.png", width = 10, height =4, units = c("in"), dpi = 300)

#Analyze specific class information
#merge two data frames
comp_data<-read.csv("Analyze_result1.csv",header=TRUE)
carnitine_data<-comp_data %>%
  filter(class=="Carnitines")
carnitine_name<-c("Carnitine C4:0","Carnitine C12:0","Carnitine C12:1","Carnitine C14:3",
                  "Carnitine C10:3","Carnitine C10:0","Carnitine C16:2","Carnitine C16:1",
                  "Carnitine C6:0","Carnitine C2:0","Carnitine","Carnitine C18:3",
                  "Carnitine C18:2","Carnitine C8:0","Carnitine C18:1","Carnitine C16:0",
                  "Carnitine C3:0","Carnitine C18:0","Carnitine C14:2","Carnitine C14:0",
                  "Carnitine C14:1","Carnitine C5:0")
carnitine_data$carnitine_name<-carnitine_name
carnitine_data<-carnitine_data %>%
  dplyr::rename(c(Low_dose=low_fc,High_dose=high_fc))
#Draw carnitine bar plot
library(tidyr)
carnitine_data_tidy<-carnitine_data %>%
  pivot_longer(names_to="type",
               values_to="fc",
               cols=c(Low_dose,High_dose))
carnitine_data_tidy$type<-factor(carnitine_data_tidy$type,levels=c('Low_dose','High_dose'))
carnitine_data_tidy$labels<-c('***','***','*','***','','*','','','','*','','',
                              '**','**','','','***','***','*','*','','',
                              '***','***','*','','','**','','*','','*','','',
                              '','','','','','','','','**','')
ggplot(carnitine_data_tidy,aes(x=carnitine_name,y=fc,fill=type))+
  scale_fill_manual(values =c("#21e1a5", "#fbe45b"))+
  scale_x_discrete(limits=c("Carnitine","Carnitine C2:0","Carnitine C3:0","Carnitine C4:0","Carnitine C5:0",
                            "Carnitine C6:0"))+
  theme(legend.position=c(0.5,0.9),legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x=NULL,y = "Fold change")+
  geom_col(position = "dodge")+
  scale_y_continuous(limits=c(0,4.5))+
  geom_hline(yintercept=1,linetype="dashed",color="#a9a9a9")+
  geom_text(aes(y=fc,label=labels),position = position_dodge(width=0.9),vjust=1.5)
ggsave("Short chain carnitines1.png", width = 7, height = 5, units = c("in"), dpi = 300)
ggplot(carnitine_data_tidy,aes(x=carnitine_name,y=fc,fill=type))+
  scale_fill_manual(values =c("#21e1a5", "#fbe45b"))+
  scale_x_discrete(limits=c("Carnitine C8:0","Carnitine C10:0","Carnitine C10:3",
                            "Carnitine C12:0","Carnitine C12:1"))+
  theme(legend.position=c(0.5,0.9),legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x=NULL,y = "Fold change")+
  geom_col(position = "dodge")+
  scale_y_continuous(limits=c(0,4.5))+
  geom_hline(yintercept=1,linetype="dashed",color="#a9a9a9")+
  geom_text(aes(y=fc,label=labels),position = position_dodge(width=0.9),vjust=1.5)
ggsave("Medium chain carnitines1.png", width = 6, height = 5, units = c("in"), dpi = 300)
ggplot(carnitine_data_tidy,aes(x=carnitine_name,y=fc,fill=type))+
  scale_fill_manual(values =c("#21e1a5", "#fbe45b"))+
  scale_x_discrete(limits=c("Carnitine C14:0","Carnitine C14:1","Carnitine C14:2",
                            "Carnitine C14:3","Carnitine C16:0","Carnitine C16:1",
                            "Carnitine C16:2","Carnitine C18:0","Carnitine C18:1",
                            "Carnitine C18:2","Carnitine C18:3"))+
  theme(legend.position=c(0.5,0.9),legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x=NULL,y = "Fold change")+
  geom_col(position = "dodge")+
  scale_y_continuous(limits=c(0,12.5))+
  geom_hline(yintercept=1,linetype="dashed",color="#a9a9a9")+
  geom_text(aes(y=fc,label=labels),position = position_dodge(width=0.9),vjust=1.5)
ggsave("Long chain carnitines1.png", width = 11, height = 5, units = c("in"), dpi = 300)

#Lipids heatmap
lipids_sum<-comp_data %>%
  select(compounds,Ctrl.1.1:PFO5DoDA.3.9,class) %>%  
  dplyr::filter(class=="Cer"|class=="CL"|class=="Co"|class=="DG"|
           class=="FFA"|class=="LPC"|class=="LPE"|class=="PC"|
           class=="PE"|class=="PI"|class=="SM"|class=="TG") %>%
  select(-class)
row.names(lipids_sum)=lipids_sum[,1]
lipids_sum[,1]<-NULL
class_info<-comp_data %>%
  select(compounds,class)
row.names(class_info)=class_info[,1]
class_info[,1]<-NULL
scaled.lipids<-lipids_sum %>%
  t %>%
  scale %>%
  t %>%
  merge(class_info,by='row.names') %>%
  group_by(class) %>%
  summarize_at(vars(Ctrl.1.1:PFO5DoDA.3.9),mean)
scaled.lipids<-as.data.frame(scaled.lipids)
row.names(scaled.lipids)=scaled.lipids[,1]
scaled.lipids[,1]<-NULL
lipid_heat_data<-as.matrix(scaled.lipids)
scaled_lipid_heat_data<-t(scale(t(lipid_heat_data)))
library(ComplexHeatmap)
library(circlize)
ha = HeatmapAnnotation(foo=anno_block(gp = gpar(fill = 2:4),
                                        labels = c("Control", "Low dose", "High dose"), 
                                        labels_gp = gpar(col = "white", fontsize = 10)))
split <- factor(c("Control", "Control","Control","Control","Control",
                  "Control","Control","Control","Control","Control",
                  "Control","Low dose","Low dose","Low dose","Low dose",
                  "Low dose","Low dose","Low dose","Low dose","Low dose",
                  "Low dose","High dose","High dose","High dose","High dose",
                  "High dose","High dose","High dose","High dose","High dose","High dose"), levels=c("Control","Low dose","High dose"))
png("lipids_heatmap1.png", width=6, height=4, units="in", res=300)
Heatmap(scaled_lipid_heat_data,cluster_columns = FALSE,
        col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Relative intensity",title_position="lefttop-rot"),
        top_annotation = ha,show_column_names = FALSE,
        column_split =split,
        cluster_column_slices = FALSE,column_title = NULL)
dev.off()
#Metabolomics heatmap
meta_data<-comp_data %>%
  dplyr::filter(platform.x=="metabolomics") %>%
  dplyr::filter(low_pvalue<0.05|high_pvalue<0.05)
meta_heat_data<-meta_data %>%
  select(Standard.name,Ctrl.1.1:PFO5DoDA.3.9)
row.names(meta_heat_data)=meta_heat_data[,1]
meta_heat_data[,1]<-NULL
m.meta_heat_data<-as.matrix(meta_heat_data)
m.meta_heat_data<-t(scale(t(m.meta_heat_data)))
lgd1 = HeatmapAnnotation(foo=anno_block(gp = gpar(fill = 2:4),
                                      labels = c("Control", "Low dose", "High dose"), 
                                      labels_gp = gpar(col = "white", fontsize = 10)))
split <- factor(c("Control", "Control","Control","Control","Control",
                  "Control","Control","Control","Control","Control",
                  "Control","Low dose","Low dose","Low dose","Low dose",
                  "Low dose","Low dose","Low dose","Low dose","Low dose",
                  "Low dose","High dose","High dose","High dose","High dose",
                  "High dose","High dose","High dose","High dose","High dose","High dose"), levels=c("Control","Low dose","High dose"))
col1=list(Class= c(FFA = "#a6cee3", Carnitines = "#1f78b4", Cholines = "#b2df8a",Steroids="#33a02c",
                'Amino acids'="#fb9a99",'Organic acids'="#e31a1c",'Fatty amides'="#fdbf6f",'Bile acids'="#ff7f00",
                Indolines="#cab2d6",Butyrophenones="#6a3d9a",LPC="#ffff99",LPE="#b15928",
                PC="#8dd3c7",PE="#ffffb3",SM="#bebada",Sphingosines="#fb8072",Sphingolipids="#80b1d3",
                Xanthines="#fdb462",Piperidinones="#b3de69",'Pyridinecarboxylic acids'="#fccde5",
                Phenylsulfates="#d9d9d9",'Pyrimidine ribonucleoside'="#bc80bd"))
lgd2=HeatmapAnnotation(Class = meta_data$class,
                            annotation_legend_param = list(
                              Class = list(
                                title = "Class",
                                at = unique(meta_data$class),
                                labels = unique(meta_data$class),
                                title_position="lefttop-rot"
                              )
                            ),show_annotation_name=FALSE,col=col1,which='row')
ht=Heatmap(m.meta_heat_data,cluster_columns = FALSE,
           col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
           heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Relative intensity",title_position="lefttop-rot"),
           top_annotation = lgd1,show_column_names = FALSE,
           column_split =split,right_annotation = lgd2,
           cluster_column_slices = FALSE,column_title = NULL)
png("metabolomics_heatmap1.png", width=12, height=20, units="in", res=300)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

#FFA structure
ffa<-comp_data %>%
  dplyr::filter(class=="FFA") %>%
  select(Standard.name,low_fc,high_fc,low_pvalue,high_pvalue) %>%
  dplyr::filter(!(Standard.name %in% c("Arachidonic acid","Trans-Vaccenic acid",
                                "2-Hydroxy-3-methylbutyric acid")))
ffa$c_num<-c(11,19,9,7,6,8,10,12,13,14,14,15,16,16,17,17,17,18,18,18,18,18,19,20,
             20,20,20,20,20,21,22,22,22,22,22,22,22,23,24,24,24,24,24,24,26)
ffa$d_num<-c(0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,
             1,2,0,1,2,3,4,0,0,1,2,3,
             4,5,0,0,1,2,3,4,5,6,0,0,
             1,2,4,5,6,0)
library(tidyr)
library(tidyverse)
ffa.low.fc<-ffa %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
ffa.low.fc<-as.matrix(ffa.low.fc)
ffa.low.log2fc<-log2(ffa.low.fc)
ffa.high.fc<-ffa %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
ffa.high.fc<-as.matrix(ffa.high.fc)
ffa.high.log2fc<-log2(ffa.high.fc)
ffa.low.p<-ffa %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
ffa.high.p<-ffa %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
cn = colnames(ffa.low.log2fc)
ffa.ht1=Heatmap(ffa.low.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
                col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                row_order = c('6', '7', '8','9', '10', '11', '12', '13', '14', '15', '16', '17', '18','19', '20', '21', '22', '23', '24', '26'),
                show_column_names = FALSE,
                bottom_annotation = HeatmapAnnotation(
                  text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                  annotation_height = max_text_width(cn)),
                heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Log2 (fold change)",title_position="lefttop-rot"),
                column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(!is.na(ffa.low.p[i, j])&ffa.low.p[i, j] <0.05&ffa.low.p[i, j] >=0.01){
                    grid.text("*", x,y)}
                  if(!is.na(ffa.low.p[i, j])&ffa.low.p[i, j] <0.01&ffa.low.p[i, j] >=0.001){
                    grid.text("**", x, y)}
                  if(!is.na(ffa.low.p[i, j])&ffa.low.p[i, j] <0.001){
                    grid.text("***", x, y)}
                })
ffa.ht2=Heatmap(ffa.high.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
                col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                row_order = c('6', '7', '8','9', '10', '11', '12', '13', '14', '15', '16', '17', '18','19', '20', '21', '22', '23', '24', '26'),
                show_column_names = FALSE,
                bottom_annotation = HeatmapAnnotation(
                  text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                  annotation_height = max_text_width(cn)),
                column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
                show_heatmap_legend = FALSE,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(!is.na(ffa.high.p[i, j])&ffa.high.p[i, j] <0.05&ffa.high.p[i, j] >=0.01){
                    grid.text("*", x,y)}
                  if(!is.na(ffa.high.p[i, j])&ffa.high.p[i, j] <0.01&ffa.high.p[i, j] >=0.001){
                    grid.text("**", x, y)}
                  if(!is.na(ffa.high.p[i, j])&ffa.high.p[i, j] <0.001){
                    grid.text("***", x, y)}
                })
ffa.ht=ffa.ht1+ffa.ht2
png("FFA_structure1.png", width=8, height=6, units="in", res=300)
draw(ffa.ht,column_title = "Association of PFO5DoDa (left: low dose; right: high dose) and FFA structure", 
     heatmap_legend_side = "right")
dev.off()

#PC structure
#When two PC molecules have the same carbon number and double bond number, we select
# the one with lower p values.
pc<-comp_data %>%
  dplyr::filter(class=="PC") %>%
  select(Standard.name,low_fc,high_fc,low_pvalue,high_pvalue) %>%
  dplyr::filter(!(Standard.name %in% c("PC O-34:2","PC O-36:4","PC O-36:5","PC O-38:5","PC O-38:6",
                                "PC(16:0e_16:0)","PC(16:0e_18:1)","PC(16:0e_20:4)","PC(16:0e_22:6)",
                                "PC(18:0e_18:2)","PC(18:0e_20:4)","PC(18:2_18:2)","PC(18:1_18:2)","PC(18:3_18:2)",
                                "PC(18:0_20:2)","PC(20:1_18:2)","PC(16:0_22:5)")))
pc$c_num<-c(30,33,34,34,35,36,32,32,34,34,
            36,36,36,38,32,35,38,35,37,39,
            36,36,38,38,40,38,40,40,37,36,
            38,40)
pc$d_num<-c(0,1,1,4,3,6,0,1,2,3,
            3,4,5,6,2,1,7,2,4,6,
            1,2,3,4,5,5,7,8,2,0,2,4)
pc.low.fc<-pc %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
pc.low.fc<-as.matrix(pc.low.fc)
pc.low.log2fc<-log2(pc.low.fc)
pc.high.fc<-pc %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
pc.high.fc<-as.matrix(pc.high.fc)
pc.high.log2fc<-log2(pc.high.fc)
pc.low.p<-pc %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
pc.high.p<-pc %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
cn = c('0','1','2','3','4','5','6','7','8')
pc.ht1=Heatmap(pc.low.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
                col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                row_order = c('30','32','33','34','35','36','37','38','39','40'),
               column_order=c('0','1','2','3','4','5','6','7','8'),
                show_column_names = FALSE,
                bottom_annotation = HeatmapAnnotation(
                  text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                  annotation_height = max_text_width(cn)),
                heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Log2 (fold change)",title_position="lefttop-rot"),
                column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(!is.na(pc.low.p[i, j])&pc.low.p[i, j] <0.05&pc.low.p[i, j] >=0.01){
                    grid.text("*", x,y)}
                  if(!is.na(pc.low.p[i, j])&pc.low.p[i, j] <0.01&pc.low.p[i, j] >=0.001){
                    grid.text("**", x, y)}
                  if(!is.na(pc.low.p[i, j])&pc.low.p[i, j] <0.001){
                    grid.text("***", x, y)}
                })
pc.ht2=Heatmap(pc.high.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
                col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                row_order = c('30','32','33','34','35','36','37','38','39','40'),
                column_order=c('0','1','2','3','4','5','6','7','8'),
                show_column_names = FALSE,
                bottom_annotation = HeatmapAnnotation(
                  text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                  annotation_height = max_text_width(cn)),
                column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
                show_heatmap_legend = FALSE,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(!is.na(pc.high.p[i, j])&pc.high.p[i, j] <0.05&pc.high.p[i, j] >=0.01){
                    grid.text("*", x,y)}
                  if(!is.na(pc.high.p[i, j])&pc.high.p[i, j] <0.01&pc.high.p[i, j] >=0.001){
                    grid.text("**", x, y)}
                  if(!is.na(pc.high.p[i, j])&pc.high.p[i, j] <0.001){
                    grid.text("***", x, y)}
                })
pc.ht=pc.ht1+pc.ht2
png("PC_structure1.png", width=8, height=6, units="in", res=300)
draw(pc.ht,column_title = "Association of PFO5DoDa (left: low dose; right: high dose) and PC structure", 
     heatmap_legend_side = "right")
dev.off()

#PE structure
pe<-comp_data %>%
  dplyr::filter(class=="PE") %>%
  select(Standard.name,low_fc,high_fc,low_pvalue,high_pvalue) %>%
  dplyr::filter(!(Standard.name %in% c("PE O-34:3","PE O-36:5","PE O-38:6","PE O-38:7",
                                "PE(16:1e_20:4)","PE(16:1e_22:4)","PE(16:1e_22:6)","PE(18:1e_18:2)",
                                "PE(18:1e_20:4)","PE(18:1e_22:4)","PE(18:1e_22:5)","PE(18:1e_22:6)",
                                "PE(18:2e_20:4)","PE(18:2e_22:6)","PE(20:1e_20:4)","PE(20:1e_22:6)",
                                "PE(16:0_18:3)","PE(18:0_22:4)","PE(18:1_20:4)")))
pe$c_num<-c(30,34,34,36,38,38,34,36,
            38,35,37,39,36,36,38,40,36,
            40,40,39,41,38,40,42)
pe$d_num<-c(0,1,2,4,5,6,3,5,
            7,2,4,6,1,2,4,6,3,
            7,8,4,6,2,4,6)
pe.low.fc<-pe %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num") 
pe.low.fc<-as.matrix(pe.low.fc)
pe.low.log2fc<-log2(pe.low.fc)
pe.high.fc<-pe %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
pe.high.fc<-as.matrix(pe.high.fc)
pe.high.log2fc<-log2(pe.high.fc)
pe.low.p<-pe %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
pe.high.p<-pe %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
cn = colnames(pe.high.fc)
pe.ht1=Heatmap(pe.low.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = c("30","34","35","36","37","38","39","40","41","42"),
               column_order=c('0','1','2','3','4','5','6','7','8'),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(cn)),
               heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Log2 (fold change)",title_position="lefttop-rot"),
               column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(pe.low.p[i, j])&pe.low.p[i, j] <0.05&pe.low.p[i, j] >=0.01){
                   grid.text("*", x,y)}
                 if(!is.na(pe.low.p[i, j])&pe.low.p[i, j] <0.01&pe.low.p[i, j] >=0.001){
                   grid.text("**", x, y)}
                 if(!is.na(pe.low.p[i, j])&pe.low.p[i, j] <0.001){
                   grid.text("***", x, y)}
               })
pe.ht2=Heatmap(pe.high.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = c("30","34","35","36","37","38","39","40","41","42"),
               column_order=c('0','1','2','3','4','5','6','7','8'),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(cn)),
               column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
               show_heatmap_legend = FALSE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(pe.high.p[i, j])&pe.high.p[i, j] <0.05&pe.high.p[i, j] >=0.01){
                   grid.text("*", x,y)}
                 if(!is.na(pe.high.p[i, j])&pe.high.p[i, j] <0.01&pe.high.p[i, j] >=0.001){
                   grid.text("**", x, y)}
                 if(!is.na(pe.high.p[i, j])&pe.high.p[i, j] <0.001){
                   grid.text("***", x, y)}
               })
(pe.ht=pe.ht1+pe.ht2)
png("PE_structure1.png", width=8, height=6, units="in", res=300)
draw(pe.ht,column_title = "Association of PFO5DoDa (left: low dose; right: high dose) and PE structure", 
     heatmap_legend_side = "right")
dev.off()

#SM structure
sm<-comp_data %>%
  dplyr::filter(class=="SM") %>%
  select(Standard.name,low_fc,high_fc,low_pvalue,high_pvalue)
sm$c_num<-c(30,32,32,33,34,35,36,36,35,
            33,34,36,39,40,40,41,41,42,42,
            42,34,43,38)
sm$d_num<-c(1,1,2,2,0,2,2,3,1,1,1,1,1,1,2,1,2,1,2,3,2,1,1)
sm.low.fc<-sm %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num") 
sm.low.fc<-as.matrix(sm.low.fc)
sm.low.log2fc<-log2(sm.low.fc)
sm.high.fc<-sm %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
sm.high.fc<-as.matrix(sm.high.fc)
sm.high.log2fc<-log2(sm.high.fc)
sm.low.p<-sm %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
sm.high.p<-sm %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
cn = colnames(sm.high.fc)
sm.ht1=Heatmap(sm.low.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = c("30","32","33", "34", "35", "36", "38","39", "40", "41", "42","43"),
               column_order=c('0','1','2','3'),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(cn)),
               heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Log2 (fold change)",title_position="lefttop-rot"),
               column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(sm.low.p[i, j])&sm.low.p[i, j] <0.05&sm.low.p[i, j] >=0.01){
                   grid.text("*", x,y)}
                 if(!is.na(sm.low.p[i, j])&sm.low.p[i, j] <0.01&sm.low.p[i, j] >=0.001){
                   grid.text("**", x, y)}
                 if(!is.na(sm.low.p[i, j])&sm.low.p[i, j] <0.001){
                   grid.text("***", x, y)}
               })
sm.ht2=Heatmap(sm.high.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = c("30","32","33", "34", "35", "36", "38","39", "40", "41", "42","43"),
               column_order=c('0','1','2','3'),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(cn)),
               column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
               show_heatmap_legend = FALSE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(sm.high.p[i, j])&sm.high.p[i, j] <0.05&sm.high.p[i, j] >=0.01){
                   grid.text("*", x,y)}
                 if(!is.na(sm.high.p[i, j])&sm.high.p[i, j] <0.01&sm.high.p[i, j] >=0.001){
                   grid.text("**", x, y)}
                 if(!is.na(sm.high.p[i, j])&sm.high.p[i, j] <0.001){
                   grid.text("***", x, y)}
               })
(sm.ht=sm.ht1+sm.ht2)
png("SM_structure1.png", width=8, height=6, units="in", res=300)
draw(sm.ht,column_title = "Association of PFO5DoDa (left: low dose; right: high dose) and SM structure", 
     heatmap_legend_side = "right")
dev.off()

#Cer structure
Cer<-comp_data %>%
  dplyr::filter(class=="Cer") %>%
  select(Standard.name,low_fc,high_fc,low_pvalue,high_pvalue) %>%
  dplyr::filter(!(Standard.name %in% c("Cer(t18:0_22:0)","Cer(t18:0_24:0)","Cer(t18:0_24:1)","Cer(d17:1_22:0)","Cer(d18:2_23:0)")))
Cer$c_num<-c(39,41,41,34,36,37,38,40,41,42,42,42,43,43,44,44,40)
Cer$d_num<-c(1,2,0,1,1,1,1,1,1,1,2,3,1,2,1,2,2)
Cer.low.fc<-Cer %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num") 
Cer.low.fc<-as.matrix(Cer.low.fc)
Cer.low.log2fc<-log2(Cer.low.fc)
Cer.high.fc<-Cer %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
Cer.high.fc<-as.matrix(Cer.high.fc)
Cer.high.log2fc<-log2(Cer.high.fc)
Cer.low.p<-Cer %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
Cer.high.p<-Cer %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
cn = colnames(Cer.high.fc)
Cer.ht1=Heatmap(Cer.low.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = c("34", "36", "37", "38", "39",  "40","41", "42", "43", "44"),
               column_order=c('0','1','2','3'),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(cn)),
               heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Log2 (fold change)",title_position="lefttop-rot"),
               column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(Cer.low.p[i, j])&Cer.low.p[i, j] <0.05&Cer.low.p[i, j] >=0.01){
                   grid.text("*", x,y)}
                 if(!is.na(Cer.low.p[i, j])&Cer.low.p[i, j] <0.01&Cer.low.p[i, j] >=0.001){
                   grid.text("**", x, y)}
                 if(!is.na(Cer.low.p[i, j])&Cer.low.p[i, j] <0.001){
                   grid.text("***", x, y)}
               })
Cer.ht2=Heatmap(Cer.high.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = c("34", "36", "37", "38", "39",  "40","41", "42", "43", "44"),
               column_order=c('0','1','2','3'),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(cn)),
               column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
               show_heatmap_legend = FALSE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(Cer.high.p[i, j])&Cer.high.p[i, j] <0.05&Cer.high.p[i, j] >=0.01){
                   grid.text("*", x,y)}
                 if(!is.na(Cer.high.p[i, j])&Cer.high.p[i, j] <0.01&Cer.high.p[i, j] >=0.001){
                   grid.text("**", x, y)}
                 if(!is.na(Cer.high.p[i, j])&Cer.high.p[i, j] <0.001){
                   grid.text("***", x, y)}
               })
(Cer.ht=Cer.ht1+Cer.ht2)
png("Cer_structure1.png", width=8, height=6, units="in", res=300)
draw(Cer.ht,column_title = "Association of PFO5DoDa (left: low dose; right: high dose) and Cer structure", 
     heatmap_legend_side = "right")
dev.off()

#DG structure
DG<-comp_data %>%
  dplyr::filter(class=="DG") %>%
  select(Standard.name,low_fc,high_fc,low_pvalue,high_pvalue) %>%
  dplyr::filter(!(Standard.name %in% c("DG(16:0_20:4)","DG(18:0_20:4)","DG(16:0_22:6)","DG(18:1_22:5)")))
DG$c_num<-c(34,34,40,36,36,38,38,40,36,38,40,36)
DG$d_num<-c(1,2,6,2,3,4,5,7,4,6,8,5)
DG.low.fc<-DG %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num") 
DG.low.fc<-as.matrix(DG.low.fc)
DG.low.log2fc<-log2(DG.low.fc)
DG.high.fc<-DG %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
DG.high.fc<-as.matrix(DG.high.fc)
DG.high.log2fc<-log2(DG.high.fc)
DG.low.p<-DG %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
DG.high.p<-DG %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
cn = colnames(DG.high.fc)
DG.ht1=Heatmap(DG.low.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
                col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                row_order = c("34","36","38","40"),
                column_order=c('1','2','3','4','5','6','7','8'),
                show_column_names = FALSE,
                bottom_annotation = HeatmapAnnotation(
                  text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                  annotation_height = max_text_width(cn)),
                heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Log2 (fold change)",title_position="lefttop-rot"),
                column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(!is.na(DG.low.p[i, j])&DG.low.p[i, j] <0.05&DG.low.p[i, j] >=0.01){
                    grid.text("*", x,y)}
                  if(!is.na(DG.low.p[i, j])&DG.low.p[i, j] <0.01&DG.low.p[i, j] >=0.001){
                    grid.text("**", x, y)}
                  if(!is.na(DG.low.p[i, j])&DG.low.p[i, j] <0.001){
                    grid.text("***", x, y)}
                })
DG.ht2=Heatmap(DG.high.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
                col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = c("34","36","38","40"),
               column_order=c('1','2','3','4','5','6','7','8'),
                show_column_names = FALSE,
                bottom_annotation = HeatmapAnnotation(
                  text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                  annotation_height = max_text_width(cn)),
                column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
                show_heatmap_legend = FALSE,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(!is.na(DG.high.p[i, j])&DG.high.p[i, j] <0.05&DG.high.p[i, j] >=0.01){
                    grid.text("*", x,y)}
                  if(!is.na(DG.high.p[i, j])&DG.high.p[i, j] <0.01&DG.high.p[i, j] >=0.001){
                    grid.text("**", x, y)}
                  if(!is.na(DG.high.p[i, j])&DG.high.p[i, j] <0.001){
                    grid.text("***", x, y)}
                })
(DG.ht=DG.ht1+DG.ht2)
png("DG_structure1.png", width=8, height=6, units="in", res=300)
draw(DG.ht,column_title = "Association of PFO5DoDa (left: low dose; right: high dose) and DG structure", 
     heatmap_legend_side = "right")
dev.off()

#PI structure
PI<-comp_data %>%
  dplyr::filter(class=="PI") %>%
  select(Standard.name,low_fc,high_fc,low_pvalue,high_pvalue) 
PI$c_num<-c(34,34,36,38,37,36,38,40,40,38,39,40)
PI$d_num<-c(1,2,4,6,4,2,4,5,6,5,4,4)
PI.low.fc<-PI %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num") 
PI.low.fc<-as.matrix(PI.low.fc)
PI.low.log2fc<-log2(PI.low.fc)
PI.high.fc<-PI %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
PI.high.fc<-as.matrix(PI.high.fc)
PI.high.log2fc<-log2(PI.high.fc)
PI.low.p<-PI %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
PI.high.p<-PI %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
cn = colnames(PI.high.fc)
PI.ht1=Heatmap(PI.low.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = c("34","36","37","38","39","40"),
               column_order=c('1','2','4','5','6'),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(cn)),
               heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Log2 (fold change)",title_position="lefttop-rot"),
               column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(PI.low.p[i, j])&PI.low.p[i, j] <0.05&PI.low.p[i, j] >=0.01){
                   grid.text("*", x,y)}
                 if(!is.na(PI.low.p[i, j])&PI.low.p[i, j] <0.01&PI.low.p[i, j] >=0.001){
                   grid.text("**", x, y)}
                 if(!is.na(PI.low.p[i, j])&PI.low.p[i, j] <0.001){
                   grid.text("***", x, y)}
               })
PI.ht2=Heatmap(PI.high.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = c("34","36","37","38","39","40"),
               column_order=c('1','2','4','5','6'),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(cn)),
               column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
               show_heatmap_legend = FALSE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(PI.high.p[i, j])&PI.high.p[i, j] <0.05&PI.high.p[i, j] >=0.01){
                   grid.text("*", x,y)}
                 if(!is.na(PI.high.p[i, j])&PI.high.p[i, j] <0.01&PI.high.p[i, j] >=0.001){
                   grid.text("**", x, y)}
                 if(!is.na(PI.high.p[i, j])&PI.high.p[i, j] <0.001){
                   grid.text("***", x, y)}
               })
(PI.ht=PI.ht1+PI.ht2)
png("PI_structure1.png", width=8, height=6, units="in", res=300)
draw(PI.ht,column_title = "Association of PFO5DoDa (left: low dose; right: high dose) and PI structure", 
     heatmap_legend_side = "right")
dev.off()

#TG structure
TG<-comp_data %>%
  dplyr::filter(class=="TG") %>%
  select(Standard.name,low_fc,high_fc,low_pvalue,high_pvalue) %>%
  dplyr::filter(!(Standard.name %in% c("TG(18:1e_16:0_16:0)","TG(18:1e_16:0_18:1)","TG(20:0e_16:0_18:2)",
                                "TG(20:0e_18:1_18:1)","TG(20:0e_18:1_18:2)","TG(16:0_18:1_18:3)",
                                "TG(18:1_18:1_18:2)","TG(18:0_20:0_20:4)","TG(14:0_18:2_22:6)",
                                "TG(16:1_18:2_20:5)","TG(16:0_18:2_22:6)","TG(18:0_18:2_22:6)",
                                "TG(18:1_18:2_22:5)","TG(14:0_18:2_20:5)","TG(16:0_14:1_22:6)",
                                "TG(16:0_18:1_22:6)","TG(16:0_14:1_18:2)","TG(18:0_18:1_18:2)",
                                "TG(18:1_18:2_18:2)","TG(19:1_16:0_20:4)","TG(16:0_12:0_22:6)",
                                "TG(16:0_18:1_22:5)","TG(18:1_18:1_20:4)","TG(16:0_18:3_22:6)",
                                "TG(16:1_18:3_22:5)","TG(16:0_20:5_22:6)","TG(16:0_22:4_22:6)",
                                "TG(16:0_22:6_22:6)","TG(18:1_22:5_22:6)")))
TG$c_num<-c(48,52,49,49,51,49,51,53,55,44,46,46,48,48,48,50,50,50,54,51,51,51,
            55,52,52,53,54,54,55,56,58,52,52,54,56,46,48,48,50,50,56,60,53,53,
            57,57,50,52,54,56,58,60,60,62,53,54,57,58,58,60,56,57,57,58,58,60,
            53,55,57,56,58,62,53,54,58,50,52,56,52,55,55,55,55,57,58,56,56,60,
            56,58,60,60,60,60,62)
TG$d_num<-c(4,8,1,2,4,3,5,7,8,2,2,1,1,2,0,1,2,3,7,1,2,3,6,2,3,1,2,5,1,1,1,4,5,
            6,6,3,3,5,4,5,9,12,2,3,6,8,0,1,4,5,7,5,10,12,4,3,2,2,8,2,7,3,4,4,9,
            3,5,7,9,8,10,14,6,8,11,6,6,10,7,2,3,4,5,7,6,2,3,6,4,5,8,9,13,7,13)
TG.low.fc<-TG %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num") 
TG.low.fc<-as.matrix(TG.low.fc)
TG.low.log2fc<-log2(TG.low.fc)
TG.high.fc<-TG %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_fc) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
TG.high.fc<-as.matrix(TG.high.fc)
TG.high.log2fc<-log2(TG.high.fc)
TG.low.p<-TG %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = low_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
TG.high.p<-TG %>%
  pivot_wider(id_cols=c_num,names_from = d_num,values_from = high_pvalue) %>%
  remove_rownames %>%
  column_to_rownames(var="c_num")
cn = colnames(TG.high.fc)
TG.ht1=Heatmap(TG.low.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = c("44","46","48","49","50","51", "52",   "53","54","55","56","57", "58", "60", "62"),
               column_order=c( "0","1","2","3", "4","5","6","7","8","9","10","11","12","13","14"),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(cn)),
               heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Log2 (fold change)",title_position="lefttop-rot"),
               column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(TG.low.p[i, j])&TG.low.p[i, j] <0.05&TG.low.p[i, j] >=0.01){
                   grid.text("*", x,y)}
                 if(!is.na(TG.low.p[i, j])&TG.low.p[i, j] <0.01&TG.low.p[i, j] >=0.001){
                   grid.text("**", x, y)}
                 if(!is.na(TG.low.p[i, j])&TG.low.p[i, j] <0.001){
                   grid.text("***", x, y)}
               })
TG.ht2=Heatmap(TG.high.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = c("44","46","48","49","50","51", "52",   "53","54","55","56","57", "58", "60", "62"),
               column_order=c( "0","1","2","3", "4","5","6","7","8","9","10","11","12","13","14"),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(cn)),
               column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
               show_heatmap_legend = FALSE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(TG.high.p[i, j])&TG.high.p[i, j] <0.05&TG.high.p[i, j] >=0.01){
                   grid.text("*", x,y)}
                 if(!is.na(TG.high.p[i, j])&TG.high.p[i, j] <0.01&TG.high.p[i, j] >=0.001){
                   grid.text("**", x, y)}
                 if(!is.na(TG.high.p[i, j])&TG.high.p[i, j] <0.001){
                   grid.text("***", x, y)}
               })
(TG.ht=TG.ht1+TG.ht2)
png("TG_structure1.png", width=8, height=6, units="in", res=300)
draw(TG.ht,column_title = "Association of PFO5DoDa (left: low dose; right: high dose) and TG structure", 
     heatmap_legend_side = "right")
dev.off()

#Completely significant metabolites and lipids heatmap
meta_data<-comp_data %>%
  dplyr::filter((low_pvalue<0.05&low_FDR<0.05)|(high_pvalue<0.05&high_FDR<0.05)) %>%
  dplyr::filter(high_fc>2|high_fc<0.5|low_fc>2|low_fc<0.5)
meta_heat_data<-meta_data %>%
  select(Standard.name,Ctrl.1.1:PFO5DoDA.3.9)
row.names(meta_heat_data)=meta_heat_data[,1]
meta_heat_data[,1]<-NULL
m.meta_heat_data<-as.matrix(meta_heat_data)
m.meta_heat_data<-t(scale(t(m.meta_heat_data)))
lgd1 = HeatmapAnnotation(foo=anno_block(gp = gpar(fill = 2:4),
                                        labels = c("Control", "Low dose", "High dose"), 
                                        labels_gp = gpar(col = "white", fontsize = 10)))
split <- factor(c("Control", "Control","Control","Control","Control",
                  "Control","Control","Control","Control","Control",
                  "Control","Low dose","Low dose","Low dose","Low dose",
                  "Low dose","Low dose","Low dose","Low dose","Low dose",
                  "Low dose","High dose","High dose","High dose","High dose",
                  "High dose","High dose","High dose","High dose","High dose","High dose"), levels=c("Control","Low dose","High dose"))
col1=list(Class= c(FFA = "#a6cee3", Carnitines = "#1f78b4", DG = "#b2df8a",Steroids="#33a02c",
                   'Amino acids'="#fb9a99",TG="#e31a1c",PI="#bc80bd",'Bile acids'="#ff7f00",
                   Cer="#cab2d6",Butyrophenones="#6a3d9a",LPC="#ffff99",LPE="#b15928",
                   PC="#8dd3c7",PE="#d9d9d9",Co="#fb8072",Sphingolipids="#80b1d3",
                   Xanthines="#fdb462",Piperidinones="#545454",'Pyridinecarboxylic acids'="#fccde5",'Fatty amides'='#8A2BE2'))
lgd2=HeatmapAnnotation(Class = meta_data$class,
                       annotation_legend_param = list(
                         Class = list(
                           title = "Class",
                           at = unique(meta_data$class),
                           labels = unique(meta_data$class),
                           title_position="lefttop-rot"
                         )
                       ),show_annotation_name=FALSE,col=col1,which='row')
ht=Heatmap(m.meta_heat_data,cluster_columns = FALSE,
           col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
           heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Relative abundance",title_position="lefttop-rot"),
           top_annotation = lgd1,show_column_names = FALSE,
           column_split =split,right_annotation = lgd2,
           cluster_column_slices = FALSE,column_title = NULL)
png("completely significant heatmap2.jpg", width=14, height=18.5, units="in", res=300)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

#Draw some boxplot in the article
ggplot(PCAdata) +
  aes(x = group, y = linolenyl.carnitine) +
  geom_boxplot() +
  theme_minimal()+
  scale_x_discrete(limits=c("Control","PFO5DoDA-low","PFO5DoDA-high"))+
  labs(y='Carnitine C18:3')+
  theme_bw()
ggsave('carntine 18_3.png',width=5, height=4, units="in", dpi=300)
ggplot(PCAdata) +
  aes(x = group, y = FA.18.3.) +
  geom_boxplot() +
  theme_minimal()+
  scale_x_discrete(limits=c("Control","PFO5DoDA-low","PFO5DoDA-high"))+
  labs(y='FA 18:3')+
  theme_bw()
ggsave('FA 18_3.png',width=5, height=4, units="in", dpi=300)

#compute carnitine ratio
PCAdata<-read.csv("PCAdata1.csv",header=TRUE,row.names = 1)
PCAdata.ratio<-PCAdata %>%
  mutate(C16.18.C0=(Palmitoylcarnitine+Stearoylcarnitine)/L.Carnitine,
         C2.C0=L.Acetylcarnitine/L.Carnitine,
         C3.C0=Propionyl.L.carnitine/L.Carnitine)
ggplot(PCAdata.ratio) +
  aes(x = group, y = C16.18.C0) +
  geom_boxplot() +
  theme_minimal()+
  scale_x_discrete(limits=c("Control","PFO5DoDA-low","PFO5DoDA-high"))+
  labs(y='C16+C18/C0')+
  theme_bw()
ggsave('C16+18_0.png',width=5, height=4, units="in", dpi=300)
ggplot(PCAdata.ratio) +
  aes(x = group, y = C2.C0) +
  geom_boxplot() +
  theme_minimal()+
  scale_x_discrete(limits=c("Control","PFO5DoDA-low","PFO5DoDA-high"))+
  labs(y='C2/C0')+
  theme_bw()
ggsave('C2_0.png',width=5, height=4, units="in", dpi=300)
ggplot(PCAdata.ratio) +
  aes(x = group, y = C3.C0) +
  geom_boxplot() +
  theme_minimal()+
  scale_x_discrete(limits=c("Control","PFO5DoDA-low","PFO5DoDA-high"))+
  labs(y='C3/C0')+
  theme_bw()
ggsave('C3_0.png',width=5, height=4, units="in", dpi=300)

#Glucocorticoid acid
Glu<-read.csv("glucocorticoid acid.csv",header=TRUE,row.names = 1)
ggplot(Glu) +
  aes(x = group, y = X11.Dehydrocorticosterone) +
  geom_boxplot() +
  theme_minimal()+
  scale_x_discrete(limits=c("Control","PFO5DoDA-low","PFO5DoDA-high"))+
  labs(y='11-Dehydrocorticosterone')+
  theme_bw()
ggplot(Glu) +
  aes(x = group, y = corticosterone) +
  geom_boxplot() +
  theme_minimal()+
  scale_x_discrete(limits=c("Control","PFO5DoDA-low","PFO5DoDA-high"))+
  labs(y='corticosterone')+
  theme_bw()+
  geom_jitter(color='blue', width=0.2,size=0.7, alpha=1)
low_glu <- Glu %>% 
  filter(group == "Control"|group == "PFO5DoDA-low")
high_glu <- Glu %>% 
  filter(group == "Control"|group == "PFO5DoDA-high")
wilcox.test(low_glu$X11.Dehydrocorticosterone ~ low_glu$group)
wilcox.test(low_glu$corticosterone ~ low_glu$group)
wilcox.test(high_glu$corticosterone ~ high_glu$group)

#Venn diagram
library(VennDiagram)
myCol <- c("#1f78b4","#b2df8a")
low_change<-comp_data %>%
  filter(low_pvalue<0.05&low_FDR<0.05) %>%
  select(Standard.name) %>%
  unlist()
high_change<-comp_data %>%
  filter(high_pvalue<0.05&high_FDR<0.05) %>%
  select(Standard.name) %>%
  unlist()
venn.diagram(
  x = list(low_change, high_change),
  category.names = c("low dose vs control","high dose vs control"),
  filename = 'venn_diagramm.png',
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  cex=0.5,
  cat.cex=0.5,
  cat.default.pos='outer',
  cat.pos = c(-25,25),
  cat.dist = c(0, 0)
)

#Lipid boxplot (we have computed 'scaled_lipid_heat_data' from lipid heatmap section)
lipid.data<-as.data.frame(scaled_lipid_heat_data) %>%
  rownames_to_column(var='category') %>%
  pivot_longer(names_to='sample',values_to = 'abundance',cols=-category)
write.csv(lipid.data,file='LipidBoxPlot.csv')
#Some modification was done on data format of csv
lipid.data<-read.csv('LipidBoxPlot.csv',header=TRUE)
compare_means(abundance~group, data = lipid.data,group.by = 'category')
my_comparisons <- list( c("low dose", "high dose"),c("control", "low dose"), c("control", "high dose"))
ggboxplot(lipid.data, x = "group", y = "abundance",
               color = "group", palette = "jco",
               add = "jitter",
               facet.by = "category", short.panel.labs = TRUE,legend='none',xlab = FALSE,
          ylab='Relative abundance')+
  stat_compare_means(label = "p.signif",method='wilcox.test',comparisons = my_comparisons,label.y = c(3, 3.6, 4.5))+
  ylim(-3.1,5.2)
ggsave('LipidBoxplot.png',width=10, height=7, units="in", dpi=300)

#Amino acids boxplot
aa.data<-comp_data %>%
  select(compounds,Ctrl.1.1:PFO5DoDA.3.9,class) %>%  
  dplyr::filter(class=="Amino acids") %>%
  select(-class) %>%
  filter(compounds %in% c('Arginine','Tryptophan','Methionine','Leucine',
                          'Tyrosine','Proline','Phenylalanine','Glutamic acid',
                          'Glutamine','Histidine'))
row.names(aa.data)=aa.data[,1]
aa.data[,1]<-NULL
aa.data<-aa.data %>%
  t %>%
  scale %>%
  t
aa.data<-as.data.frame(aa.data) %>%
  rownames_to_column(var='compound') %>%
  pivot_longer(names_to='sample',values_to = 'abundance',cols=-compound)
write.csv(aa.data,file='AminoAcidBoxPlot.csv')
#Some modification was done on data format of csv
aa.data<-read.csv('AminoAcidBoxPlot.csv',header=TRUE)
aa.data$compound <- factor(aa.data$compound, c('Arginine','Tryptophan','Methionine','Leucine',
                                                 'Tyrosine','Proline','Phenylalanine','Glutamic acid',
                                                 'Glutamine','Histidine'))
compare_means(abundance~group, data = aa.data,group.by = 'compound')
my_comparisons <- list( c("low dose", "high dose"),c("control", "low dose"), c("control", "high dose"))
ggboxplot(aa.data, x = "group", y = "abundance",
          color = "group", palette = "jco",
          add = "jitter",
          facet.by = "compound", short.panel.labs = TRUE,legend='none',xlab = FALSE,
          ylab='Relative abundance',ncol=5,cex.lab=1)+
  stat_compare_means(label = "p.signif",method='wilcox.test',comparisons = my_comparisons,label.y = c(3.5, 4.1, 5))+
  ylim(-2,6)+
  theme(axis.text.x = element_text(size = 8))
ggsave('AABoxplot.png',width=8.6, height=3.8, units="in", dpi=300)

#Coenzyme plot
co_data<-comp_data %>%
  filter(class=="Co")
co_data<-co_data %>%
  rename(c('Low dose'=low_fc,'High dose'=high_fc))
co_name<-c('Coenzyme Q10','Coenzyme Q8','Coenzyme Q9')
co_data$co_name<-co_name
co_data_tidy<-co_data %>%
  pivot_longer(names_to="type",
               values_to="fc",
               cols=c('Low dose','High dose'))
co_data_tidy$type<-factor(co_data_tidy$type,levels=c('Low dose','High dose'))
co_data_tidy$labels<-c('****','****','***','***','****','****')
ggplot(co_data_tidy,aes(x=co_name,y=fc,fill=type))+
  scale_fill_manual(values =c("#21e1a5", "#fbe45b"))+
  scale_x_discrete(limits=c('Coenzyme Q8','Coenzyme Q9','Coenzyme Q10'))+
  theme(legend.position=c(0.5,0.9),legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x=NULL,y = "Fold change")+
  geom_col(position = "dodge")+
  scale_y_continuous(limits=c(0,2.5))+
  geom_hline(yintercept=1,linetype="dashed",color="#a9a9a9")+
  geom_text(aes(y=fc,label=labels),position = position_dodge(width=0.9),vjust=1.5)
ggsave("Coenzyme.png", width = 4, height = 2.5, units = c("in"), dpi = 300)



# 2023.5.24 要求分析嘌呤，提峰后对数据单独分析
library(ggplot2)
library(dplyr)
library(ggsignif)

# Read the CSV file
data <- read.csv("20230524purine append.csv")

# Group the compounds
data <- data %>%
  mutate(group = factor(group))

# Create the plot
plot <- ggplot(data, aes(x = group, y = guanine, fill = group)) +
  geom_boxplot() +
  labs(x = "Group", y = "Guanine") +
  ggtitle("Compound Distribution") +
  theme_bw()

# Add significance indication
plot <- plot +
  geom_signif(comparisons = list(c("Control", "PFO5DoDA 3"), c("Control", "PFO5DoDA 2"), c("Control", "PFO4DA 2"), c("Control", "PFO4DA 3")),
              map_signif_level =c("***"=0.001, "**"=0.01, "*"=0.05),
              textsize = 4,
              vjust = -0.5,
              test = "wilcox.test",
              step_increase = 0.05,
              tip_length = 0)

# Display the plot
plot
ggsave("20230523/guanine.png", plot, width = 8, height = 6, dpi = 300)


