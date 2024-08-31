library(tidyverse)
library(factoextra)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)    # To add significance in ggplot
library(VennDiagram)

# Install fonts
install.packages("extrafont")
library(extrafont)
font_import()
loadfonts(device="win")       #Register fonts for Windows bitmap output
fonts()                       #vector of font family names

# Another way to choose font
windowsFonts()
names(wf[wf=="TT Times New Roman"])

# Fig. S2A
# QC violin plot
metabolomics_qc<-read.csv("metabolomics_violin_plot.csv",header =TRUE)
#par(mfrow = c(1, 2))
pdf('figures2a1.pdf',width=3.6,height=3.4)
ggplot(metabolomics_qc,aes(x=mode,y=cv,fill=mode))+
  geom_violin(adjust=0.7)+
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize=0.7,binwidth = 0.4,show.legend=FALSE,fill="black")+
  scale_fill_manual(values=c('gray','red'))+
  labs(title="Metabolomics", x=NULL,y = "Coefficient of variation(%)")+
  theme(legend.position="none",axis.text=element_text(size=10,face='bold'))+
  scale_x_discrete(limits=c("pos", "neg"))+
  theme(text = element_text(family = "serif"))
dev.off()
#ggsave("FigureS2A_MetabolomicsQC.png", width = 3.6, height =3.4, units = c("in"), dpi = 300)
lipidomics_qc<-read.csv("lipidomics_violin_plot.csv",header =TRUE)
pdf('figures2a2.pdf',width=3.6,height=3.4)
ggplot(lipidomics_qc,aes(x=mode,y=cv,fill=mode))+
  geom_violin(adjust=0.7)+
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize=0.7,binwidth = 0.4,show.legend=FALSE,fill="black")+
  scale_fill_manual(values=c('gray','red'))+
  labs(title="Lipidomics", x=NULL,y = "Coefficient of variation(%)")+
  theme(legend.position="none", text = element_text(family = "serif"),axis.text=element_text(size=10,face='bold'))+
  scale_x_discrete(limits=c("pos", "neg"))
dev.off()
#ggsave("FigureS2A_LipidomicsQC.png", width = 3.6, height =3.4, units = c("in"), dpi = 300)

# Fig. S2B
compounds_sum<-read.csv("compounds_sum1.csv",header=TRUE)
x_p<-c("Carnitines","Amino acids","Bile acids","Organic acids","Others","LPC","LPE","PC","PE","SM","FFA","Cer","CL","CoQ","DG","PI","TG")
y_p<-c(22,16,8,7,20,40,18,49,43,23,48,22,8,3,16,12,124)
pdf('figures2b.pdf',width=5,height=3.5)
ggplot(compounds_sum,aes(x=class,y=count,fill=platform))+
  geom_bar(stat="identity")+
  scale_fill_manual(values =c('red','blue'))+
  scale_x_discrete(limits=c("Carnitines","Amino acids","Bile acids","Organic acids","Others","LPC","LPE","PC","PE","SM","FFA","Cer","CL","CoQ","DG","PI","TG"))+
  theme(legend.position=c(0.5,0.9),legend.title = element_blank(),legend.text = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "serif"),axis.text=element_text(size=10,face='bold'),
        axis.title=element_text(size=13,face="bold"),axis.text.x = element_text(angle = 45,hjust=1))+
  labs(x=NULL,y = "Count")+
  annotate(geom="text", x=x_p, y=y_p+3.8, label=y_p,family='serif',size=4.5)
dev.off()
#ggsave("FigureS2B_CompoundClass1.png", width = 5, height = 3.5, units = c("in"), dpi = 300)

# Fig.2A
# PCA score plot
PCAdata<-read.csv("PCAdata1.csv",header=TRUE,row.names = 1)
PCAdata.active<-PCAdata[,2:480]
res.pca <- prcomp(PCAdata.active, scale = TRUE,)
fviz_eig(res.pca)
PCAdata$group<-factor(PCAdata$group,levels=c('Control','PFO5DoDA-low','PFO5DoDA-high'), labels = c("Ctrl","2 μg/kg/d","10 μg/kg/d"))
group <- as.factor(PCAdata$group)
pdf('figure2a_1.pdf',width=6,height=5)
fviz_pca_ind(res.pca,
             col.ind =group, # color by group
             geom ="point",legend.title="Group",
             palette=c("green", "blue", "red"), addEllipses=TRUE,ellipse.level=0.6,
             invisible="quali",pointsize=2)+
  theme(text = element_text(family = "serif"),axis.title=element_text(size=18),legend.title =element_text(size=18,face="bold"),
        legend.text = element_text(size=17))+
  labs(title=(''))
dev.off()
ggsave("Figure2A.png", bg='white',width = 6, height = 5, units = c("in"), dpi = 300)

#volcano plot
compound_data<-read.csv("Analyze_result1.csv",header=TRUE)
compound_data <- compound_data %>% 
  mutate(
    log2low_fc=log2(low_fc),
    log2high_fc=log2(high_fc)
  )
# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  compound_data$log2low_fc < -0.5, 'blue',
  ifelse(compound_data$log2low_fc > 0.5, 'red',
         'gray39'))
keyvals[is.na(keyvals)] <- 'gray39'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'gray39'] <- 'mid'
names(keyvals)[keyvals == 'blue'] <- 'low'
pdf('figures2c1.pdf',width=4.6,height=3.7)
EnhancedVolcano(compound_data,
                    lab =compound_data$Standard.name, 
                    x = 'log2low_fc',
                    y = 'low_pvalue',
                    axisLabSize = 10,
                    title = 'Group: 2 μg/kg/d vs Ctrl',
                    titleLabSize = 14,
                    pCutoff = 0.05,
                    FCcutoff = 0.5,
                    pointSize = 0.75,
                    labSize = 3,
                    colCustom = keyvals,
                    colAlpha = 0.65,
                    legendPosition = 'none',
                    subtitle =NULL,
                    caption = bquote(~Log[2]~ "fold change cutoff, 0.5; p-value cutoff, 0.05"),
                    captionLabSize = 11,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
)+
  theme(text = element_text(family = "serif"))+
  ggplot2::coord_cartesian(ylim=c(0, 6),xlim=c(-5.5,3))
dev.off()
#ggsave("FigureS2C_1.png", width = 4.6, height =3.7, units = c("in"), dpi = 300)
keyvals_h <- ifelse(
  compound_data$log2high_fc < -0.5, 'blue',
  ifelse(compound_data$log2high_fc > 0.5, 'red',
         'gray39'))
keyvals_h[is.na(keyvals_h)] <- 'gray39'
names(keyvals_h)[keyvals_h == 'red'] <- 'high'
names(keyvals_h)[keyvals_h == 'gray39'] <- 'mid'
names(keyvals_h)[keyvals_h == 'blue'] <- 'low'
pdf('figures2c2.pdf',width=4.6,height=3.7)
EnhancedVolcano(compound_data,
                   lab =compound_data$Standard.name, 
                   x = 'log2high_fc',
                   y = 'high_pvalue',
                   axisLabSize = 10,
                   title = 'Group: 10 μg/kg/d vs Ctrl',
                   titleLabSize = 14,
                   pCutoff = 0.05,
                   FCcutoff = 0.5,
                   pointSize = 0.75,
                   labSize = 3,
                colCustom = keyvals_h,
                colAlpha = 0.65,
                   legendPosition = 'none',
                   subtitle =NULL,
                   caption = bquote(~Log[2]~ "fold change cutoff, 0.5; p-value cutoff, 0.05"),
                   captionLabSize = 11,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
)+
  theme(text = element_text(family = "serif"))+
  ggplot2::coord_cartesian(ylim=c(0, 6),xlim=c(-5.3,4.3))
dev.off()
#ggsave("FigureS2C_2.png", width = 4.6, height =3.7, units = c("in"), dpi = 300)

#FFA structure
comp_data<-read.csv("Analyze_result1.csv",header=TRUE)
ffa<-comp_data %>%
  filter(class=="FFA") %>%
  select(Standard.name,low_fc,high_fc,low_pvalue,high_pvalue) %>%
  filter(!(Standard.name %in% c("Arachidonic acid","Trans-Vaccenic acid",
                                       "2-Hydroxy-3-methylbutyric acid")))
ffa$c_num<-c(11,19,9,7,6,8,10,12,13,14,14,15,16,16,17,17,17,18,18,18,18,18,19,20,
             20,20,20,20,20,21,22,22,22,22,22,22,22,23,24,24,24,24,24,24,26)
ffa$d_num<-c(0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,
             1,2,0,1,2,3,4,0,0,1,2,3,
             4,5,0,0,1,2,3,4,5,6,0,0,
             1,2,4,5,6,0)
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
a=c(26,24:6)
ffa.ht1=Heatmap(ffa.low.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
                col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                row_order = as.character(a),
                show_column_names = FALSE,
                bottom_annotation = HeatmapAnnotation(
                  text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right",gp=gpar(fontfamily = "serif")),
                  annotation_height = max_text_width(cn)),
                heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Log2 (fold change)",title_position="lefttop-rot",labels_gp=gpar(fontfamily = "serif",fontsize=13),title_gp=gpar(fontface='bold',fontfamily = "serif",fontsize=13)),
                column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(!is.na(ffa.low.p[i, j])&ffa.low.p[i, j] <0.05&ffa.low.p[i, j] >=0.01){
                    grid.text("*", x,y)}
                  if(!is.na(ffa.low.p[i, j])&ffa.low.p[i, j] <0.01&ffa.low.p[i, j] >=0.001){
                    grid.text("**", x, y)}
                  if(!is.na(ffa.low.p[i, j])&ffa.low.p[i, j] <0.001){
                    grid.text("***", x, y)}
                },row_title_gp=gpar(fontfamily = "serif",fontface='bold',fontsize=16),column_title_gp = gpar(fontfamily = "serif",fontface='bold',fontsize=16),
                row_names_gp = gpar(fontfamily = "serif",fontsize=15))
ffa.ht2=Heatmap(ffa.high.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
                col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                row_order = as.character(a),
                show_column_names = FALSE,
                bottom_annotation = HeatmapAnnotation(
                  text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right",gp=gpar(fontfamily = "serif",fontsize=13)),
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
                },row_title_gp=gpar(fontfamily = "serif",fontsize=16),column_title_gp = gpar(fontfamily = "serif",fontface='bold',fontsize=16),
                row_names_gp = gpar(fontfamily = "serif",fontsize=15))
(ffa.ht=ffa.ht1+ffa.ht2)
#png("Figure2C.png", width=7.2, height=6, units="in", res=300)
pdf('figure2e_1.pdf',width=7.2,height=6)
draw(ffa.ht,column_title = "Association of PFO5DoDA (left: 2 μg/kg/d; right: 10 μg/kg/d) and FA structure", 
     heatmap_legend_side = "right",
     column_title_gp = gpar(fontfamily = "serif",fontface='bold',fontsize=13.6))
dev.off()

#TG structure
TG<-comp_data %>%
  filter(class=="TG") %>%
  select(Standard.name,low_fc,high_fc,low_pvalue,high_pvalue) %>%
  filter(!(Standard.name %in% c("TG(18:1e_16:0_16:0)","TG(18:1e_16:0_18:1)","TG(20:0e_16:0_18:2)",
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
a=c(62,60,58:48,46,44)
TG.ht1=Heatmap(TG.low.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = as.character(a),
               column_order=c( "0","1","2","3", "4","5","6","7","8","9","10","11","12","13","14"),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right",gp=gpar(fontfamily = "serif")),
                 annotation_height = max_text_width(cn)),
               heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Log2 (fold change)",title_position="lefttop-rot",labels_gp=gpar(fontfamily = "serif",fontsize=13.3),title_gp=gpar(fontface='bold',fontfamily = "serif",fontsize=13.3)),
               column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(TG.low.p[i, j])&TG.low.p[i, j] <0.05&TG.low.p[i, j] >=0.01){
                   grid.text("*", x,y)}
                 if(!is.na(TG.low.p[i, j])&TG.low.p[i, j] <0.01&TG.low.p[i, j] >=0.001){
                   grid.text("**", x, y)}
                 if(!is.na(TG.low.p[i, j])&TG.low.p[i, j] <0.001){
                   grid.text("***", x, y)}
               },row_title_gp=gpar(fontfamily = "serif",fontface='bold',fontsize=16),column_title_gp = gpar(fontfamily = "serif",fontface='bold',fontsize=16),
               row_names_gp = gpar(fontfamily = "serif",fontsize=15))
TG.ht2=Heatmap(TG.high.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = as.character(a),
               column_order=c( "0","1","2","3", "4","5","6","7","8","9","10","11","12","13","14"),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right",gp=gpar(fontfamily = "serif")),
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
               },row_title_gp=gpar(fontfamily = "serif",fontface='bold',fontsize=16),column_title_gp = gpar(fontfamily = "serif",fontface='bold',fontsize=16),
               row_names_gp = gpar(fontfamily = "serif",fontsize=15))
(TG.ht=TG.ht1+TG.ht2)
#png("Figure2D.png", width=8, height=6.67, units="in", res=300)
pdf('figures3a.pdf',width=8,height=6.67)
draw(TG.ht,column_title = "Association of PFO5DoDA (left: 2 μg/kg/d; right: 10 μg/kg/d) and TG structure", 
     heatmap_legend_side = "right",
     column_title_gp = gpar(fontfamily = "serif",fontface='bold',fontsize=13.8))
dev.off()

#DG structure
comp_data<-read.csv("Analyze_result1.csv",header=TRUE)
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
a=c(40,38,36,34)
DG.ht1=Heatmap(DG.low.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order = as.character(a),
               column_order=c('1','2','3','4','5','6','7','8'),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right",gp=gpar(fontfamily = "serif")),
                 annotation_height = max_text_width(cn)),
               heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Log2 (fold change)",title_position="lefttop-rot",labels_gp=gpar(fontfamily = "serif",fontsize=13.3),title_gp=gpar(fontface='bold',fontfamily = "serif",fontsize=13.3)),
               column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(DG.low.p[i, j])&DG.low.p[i, j] <0.05&DG.low.p[i, j] >=0.01){
                   grid.text("*", x,y)}
                 if(!is.na(DG.low.p[i, j])&DG.low.p[i, j] <0.01&DG.low.p[i, j] >=0.001){
                   grid.text("**", x, y)}
                 if(!is.na(DG.low.p[i, j])&DG.low.p[i, j] <0.001){
                   grid.text("***", x, y)}
               },row_title_gp=gpar(fontfamily = "serif",fontface='bold',fontsize=16),column_title_gp = gpar(fontfamily = "serif",fontface='bold',fontsize=16),
               row_names_gp = gpar(fontfamily = "serif",fontsize=15))
DG.ht2=Heatmap(DG.high.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               row_order =as.character(a),
               column_order=c('1','2','3','4','5','6','7','8'),
               show_column_names = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right",gp=gpar(fontfamily = "serif")),
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
               },row_title_gp=gpar(fontfamily = "serif",fontface='bold',fontsize=16),column_title_gp = gpar(fontfamily = "serif",fontface='bold',fontsize=16),
               row_names_gp = gpar(fontfamily = "serif",fontsize=15))
(DG.ht=DG.ht1+DG.ht2)
pdf('figures3b.pdf',width=8,height=6.67)
#png("FigureS2E.png", width=8, height=6.67, units="in", res=300)
draw(DG.ht,column_title = "Association of PFO5DoDA (left: 2 μg/kg/d; right: 10 μg/kg/d) and DG structure", 
     heatmap_legend_side = "right",
     column_title_gp = gpar(fontfamily = "serif",fontface='bold',fontsize=13.8))
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
a=c(44:36,34)
Cer.ht1=Heatmap(Cer.low.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
                col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                row_order = as.character(a),
                column_order=c('0','1','2','3'),
                show_column_names = FALSE,
                bottom_annotation = HeatmapAnnotation(
                  text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right",gp=gpar(fontfamily = "serif")),
                  annotation_height = max_text_width(cn)),
                heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Log2 (fold change)",title_position="lefttop-rot",labels_gp=gpar(fontfamily = "serif",fontsize=13.3),title_gp=gpar(fontface='bold',fontfamily = "serif",fontsize=13.3)),
                column_title = "Number of Double bonds", column_title_side = "bottom",row_title = "Number of Carbon atoms",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(!is.na(Cer.low.p[i, j])&Cer.low.p[i, j] <0.05&Cer.low.p[i, j] >=0.01){
                    grid.text("*", x,y)}
                  if(!is.na(Cer.low.p[i, j])&Cer.low.p[i, j] <0.01&Cer.low.p[i, j] >=0.001){
                    grid.text("**", x, y)}
                  if(!is.na(Cer.low.p[i, j])&Cer.low.p[i, j] <0.001){
                    grid.text("***", x, y)}
                },row_title_gp=gpar(fontfamily = "serif",fontface='bold',fontsize=16),column_title_gp = gpar(fontfamily = "serif",fontface='bold',fontsize=16),
                row_names_gp = gpar(fontfamily = "serif",fontsize=15))
Cer.ht2=Heatmap(Cer.high.log2fc,cluster_columns = FALSE,cluster_rows = FALSE,
                col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                row_order = as.character(a),
                column_order=c('0','1','2','3'),
                show_column_names = FALSE,
                bottom_annotation = HeatmapAnnotation(
                  text = anno_text(cn, rot = 0, offset = unit(1, "npc"), just = "right",gp=gpar(fontfamily = "serif")),
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
                },row_title_gp=gpar(fontfamily = "serif",fontface='bold',fontsize=16),column_title_gp = gpar(fontfamily = "serif",fontface='bold',fontsize=16),
                row_names_gp = gpar(fontfamily = "serif",fontsize=15))
(Cer.ht=Cer.ht1+Cer.ht2)
pdf('figures3c.pdf',width=8,height=6)
#png("FigureS2f.png", width=8, height=6, units="in", res=300)
draw(Cer.ht,column_title = "Association of PFO5DoDA (left: 2 μg/kg/d; right: 10 μg/kg/d) and Cer structure", 
     heatmap_legend_side = "right",
     column_title_gp = gpar(fontfamily = "serif",fontface='bold',fontsize=13.8))
dev.off()


#Completely significant metabolites and lipids heatmap
meta_data<-comp_data %>%
  filter((low_pvalue<0.05&low_FDR<0.05)|(high_pvalue<0.05&high_FDR<0.05)) %>%
  filter(high_fc>2|high_fc<0.5|low_fc>2|low_fc<0.5)
#write.csv(meta_data,'RowOrder.csv',row.names = FALSE)  # Reoder in Excel
meta_data<-read.csv('RowOrder.csv',header = TRUE)
meta_heat_data<-meta_data %>%
  select(Standard.name,Ctrl.1.1:PFO5DoDA.3.9)
row.names(meta_heat_data)=meta_heat_data[,1]
meta_heat_data[,1]<-NULL
m.meta_heat_data<-as.matrix(meta_heat_data)
m.meta_heat_data<-t(scale(t(m.meta_heat_data)))
lgd1 = HeatmapAnnotation(foo=anno_block(gp = gpar(fill = 2:4),
                                        labels = c("Ctrl", "2 μg/kg/d", "10 μg/kg/d"), 
                                        labels_gp = gpar(col = "black", fontsize = 18,fontfamily='serif',fontface='bold')))
split <- factor(c("Control", "Control","Control","Control","Control",
                  "Control","Control","Control","Control","Control",
                  "Control","Low dose","Low dose","Low dose","Low dose",
                  "Low dose","Low dose","Low dose","Low dose","Low dose",
                  "Low dose","High dose","High dose","High dose","High dose",
                  "High dose","High dose","High dose","High dose","High dose","High dose"), levels=c("Control","Low dose","High dose"))
col1=list(Class= c(FA = "#276419", Carnitines = "#1f78b4", DG = "#b2df8a",Steroids="#33a02c",
                   'Amino acids'="#fb9a99",PI="#fdbf6f",'Bile acids'="#ff7f00",
                   Cer="#cab2d6",Butyrophenones="#6a3d9a",LPC="#ffff99",LPE="#b15928",
                   PC="#8dd3c7",PE="#bebada",CoQ="#fb8072",TG="#2166ac",
                   Xanthines="#fdb462",Piperidinones="#545454",'Pyridinecarboxylic acids'="#fccde5",'Fatty amides'='#8e0152'))
lgd2=HeatmapAnnotation(Class = meta_data$class,
                       annotation_legend_param = list(
                         Class = list(
                           title = "Class",
                           at = unique(meta_data$class),
                           labels = unique(meta_data$class),
                           labels_gp=gpar(fontfamily = "serif",fontsize=16),title_gp=gpar(fontface='bold',fontfamily = "serif",fontsize=18),
                           legend_gp= gpar(fill = 1:10), nrow=4
                         )
                       ),show_annotation_name=FALSE,col=col1,which='row')
ht=Heatmap(m.meta_heat_data,cluster_columns = FALSE,cluster_rows = FALSE,
           col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
           heatmap_legend_param =list(at=c(-2,-1,0,1,2),title="Relative abundance",direction = "horizontal",
                                      labels_gp=gpar(fontfamily = "serif",fontsize=16),title_gp=gpar(fontface='bold',fontfamily = "serif",fontsize=18)),
           top_annotation = lgd1,show_column_names = FALSE,
           column_split =split,right_annotation = lgd2,
           cluster_column_slices = FALSE,column_title = NULL,row_names_gp = gpar(fontfamily = "serif"))
ht
#png("FigureS2D.jpg", width=8, height=22, units="in", res=300)
pdf('figure2c.pdf',width=8,height=22)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()

#compute C2/CO ratio
PCAdata<-read.csv("PCAdata1.csv",header=TRUE,row.names = 1)
PCAdata.ratio<-PCAdata %>%
  mutate(C2.C0=L.Acetylcarnitine/L.Carnitine)
PCAdata.ratio$group<-factor(PCAdata.ratio$group,levels=c('Control','PFO5DoDA-low','PFO5DoDA-high'), 
                            labels = c("Ctrl","2 μg/kg/d","10 μg/kg/d"))
compare_means(C2.C0~group, data = PCAdata.ratio)
my_comparisons <- list( c("2 μg/kg/d", "10 μg/kg/d"),c("Ctrl", "2 μg/kg/d"), c("Ctrl", "10 μg/kg/d"))
pdf('figures3d.pdf',width=1.8,height=2.8)
ggplot(PCAdata.ratio,aes(x = group, y = C2.C0,fill=group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x=NULL,y='C2/C0')+
  #geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual('Group',values=c('white','red','gray'))+
  theme_classic(base_family = 'serif')+
  theme(legend.position="none",axis.text.x = element_text(angle = 45,hjust=1))+
  stat_compare_means(label = "p.signif",method='wilcox.test',comparisons = my_comparisons,label.y = c(1.3, 1.5, 1.7))
dev.off()
#ggsave('FigureS2J.png',width=3, height=2.8, units="in", dpi=300)

#Venn diagram
comp_data<-read.csv("Analyze_result1.csv",header=TRUE)
myCol <- c("red","blue")
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
  category.names = c("2 μg/kg/d vs Ctrl","10 μg/kg/d vs Ctrl"),
  filename = 'Figure2D.tiff',
  output=TRUE,
  # Output features
  imagetype="tiff" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  cex=0.8,
  cat.cex=0.65,
  cat.default.pos='outer',
  cat.pos = c(-25,25),
  cat.dist = c(0, 0)
)

# Lipid boxplot (we have computed 'scaled_lipid_heat_data' from lipid heatmap section)
lipid.data<-read.csv('LipidBoxPlot.csv',header=TRUE)
lipid.data$group<-factor(lipid.data$group,levels=c('control','low dose','high dose'), 
                            labels = c("Ctrl","2 μg/kg/d","10 μg/kg/d"))
compare_means(abundance~group, data = lipid.data,group.by = 'category')
my_comparisons <- list( c("2 μg/kg/d", "10 μg/kg/d"),c("Ctrl", "2 μg/kg/d"), c("Ctrl", "10 μg/kg/d"))
pdf('figures2d.pdf',width=10,height=7)
ggboxplot(lipid.data, x = "group", y = "abundance",
           fill = 'group',
          #add = "jitter",
          facet.by = "category", short.panel.labs = TRUE,legend='none',xlab = FALSE,
          ylab='Relative abundance',panel.labs.font = list(face='bold'))+
  stat_compare_means(label = "p.signif",method='wilcox.test',comparisons = my_comparisons,
                     label.y = c(3, 3.6, 4.5),
                     symnum.args = list(cutpoints = c(0,  0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")))+
  ylim(-3.1,5.2)+
  scale_fill_manual(values=c('white','red','gray'))+
  theme(text = element_text(family = "serif"))
dev.off()
#ggsave('FigureS2F.png',width=10, height=7, units="in", dpi=300)

# Amino acids boxplot
aa.data<-read.csv('AminoAcidBoxPlot.csv',header=TRUE)
aa.data$compound <- factor(aa.data$compound, c('Arginine','Tryptophan','Methionine','Leucine',
                                               'Tyrosine','Proline','Phenylalanine','Glutamic acid',
                                               'Glutamine','Histidine'))
aa.data$group<-factor(aa.data$group,levels=c('control','low dose','high dose'), 
                         labels = c("Ctrl","2 μg/kg/d","10 μg/kg/d"))
compare_means(abundance~group, data = aa.data,group.by = 'compound')
my_comparisons <- list( c("2 μg/kg/d", "10 μg/kg/d"),c("Ctrl", "2 μg/kg/d"), c("Ctrl", "10 μg/kg/d"))
pdf('figure2d_aabox.pdf',width=6,height=10)
ggboxplot(aa.data, x = "group", y = "abundance",
          fill = 'group',panel.labs.font = list(face='bold'),
          add = "none", add.params=list(color='gray18',size=0.5),
          facet.by = "compound", short.panel.labs = TRUE,legend='none',xlab = FALSE,
          ylab='Relative abundance',ncol=2,cex.lab=1)+
  stat_compare_means(label = "p.signif",method='wilcox.test',comparisons = my_comparisons,label.y = c(3.5, 4.1, 5),
                     symnum.args = list(cutpoints = c(0,  0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")))+
  ylim(-2,6)+
  theme(axis.text.x = element_text(size = 8),text = element_text(family = "serif"))+
  scale_fill_manual(values=c('white','red','gray'))
dev.off()
#ggsave('Figure2D(aa boxplot)1.png',width=6, height=8, units="in", dpi=300)

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
co_data_tidy$type<-factor(co_data_tidy$type,levels=c('Low dose','High dose'),labels = c("2 μg/kg/d","10 μg/kg/d"))
co_data_tidy$labels<-c('***','***','***','***','***','***')
pdf('figures3e.pdf',width=2,height=3)
ggplot(co_data_tidy,aes(x=co_name,y=fc,fill=type))+
  scale_fill_manual(values =c("gray", "red"))+
  scale_x_discrete(limits=c('Coenzyme Q8','Coenzyme Q9','Coenzyme Q10'))+
  theme(legend.position=c(0.1,0.9),legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "serif"),axis.text.x = element_text(angle = 45,hjust=1))+
  labs(x=NULL,y = "Fold change")+
  geom_col(position = "dodge")+
  scale_y_continuous(limits=c(0,2.5))+
  geom_hline(yintercept=1,linetype="dashed",color="gray32")+
  geom_text(aes(y=fc,label=labels),position = position_dodge(width=0.9),vjust=1.5)
dev.off()
#ggsave("FigureS2H.png", width = 3.5, height = 2.5, units = c("in"), dpi = 300)

# Figure5A
#volcano plot
gene_data<-read.csv("./FHY/Genes.csv",header=TRUE)
# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  gene_data$Log2FC < -2, 'blue',
  ifelse(gene_data$Log2FC > 1.99999, 'red',
         'gray39'))
keyvals[is.na(keyvals)] <- 'gray39'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'gray49'] <- 'mid'
names(keyvals)[keyvals == 'blue'] <- 'low'
pdf('figure4a.pdf',width=4.6,height=3.7)
EnhancedVolcano(gene_data,
                lab =gene_data$GeneName, 
                x = 'Log2FC',
                y = 'Pvalue',
                title = '',
                selectLab = c('Lpin1','Ppp1r3b','Tcim','Fkbp5',
                              'Mthfr','Tsc22d3','Tat','Slc45a3',
                              'Arrdc2','Slco1a4','Igfbp1','Zbtb16',
                              'Srebf1','Pctp','Krt23','Osbpl3','Cd36','Clstn3β'),
                axisLabSize = 10,
                titleLabSize = 12,
                pCutoff = 10e-8,
                FCcutoff = 2,
                pointSize = 0.6,
                labSize = 3,
                colCustom = keyvals,
                colAlpha = 0.65,
                legendPosition = 'none',
                subtitle =NULL,
                caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, "~10^-8),
                captionLabSize = 9,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = TRUE
)+
  theme(text = element_text(family = "serif"))+
  ggplot2::coord_cartesian(xlim=c(-7,7))
dev.off()
ggsave("Figure4A.png", width = 4.6, height =3.7, units = c("in"), dpi = 300)

# Figure 4B
F4b_data<-read.csv("./FHY/Figure4B.csv",header = TRUE)
row.names(F4b_data)=F4b_data[,1]
F4b_data[,1]<-NULL
m.F4b_data<-as.matrix(F4b_data)
lgd1 = HeatmapAnnotation(foo=anno_block(gp = gpar(fill = 2:4),
                                        labels = c("Ctrl","10 μg/kg/d"), 
                                        labels_gp = gpar(col = "black", fontsize = 18,fontfamily='serif',fontface='bold')))
split <- factor(c("Control", "Control","Control","High dose","High dose","High dose"), levels=c("Control","High dose"))
ht=Heatmap(m.F4b_data,cluster_columns = FALSE,cluster_rows = TRUE,
           col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
           heatmap_legend_param =list(at=c(-1.5,0,1.5),title="Relative abundance",direction = "horizontal",
                                      labels_gp=gpar(fontfamily = "serif",fontsize=16),title_gp=gpar(fontface='bold',fontfamily = "serif",fontsize=18)),
           top_annotation = lgd1,show_column_names = FALSE,
           column_split =split,
           cluster_column_slices = FALSE,column_title = NULL,row_names_gp = gpar(fontfamily = "serif",fontsize=18))
ht
#png("Figure4B.jpg", width=5, height=5, units="in", res=300)
pdf('figure4b.pdf',width=5,height=5)
draw(ht,  heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()

# Figure 4C
F4c_data<-read.csv("./FHY/Figure4C.csv",header = TRUE)
row.names(F4c_data)=F4c_data[,1]
F4c_data[,1]<-NULL
m.F4c_data<-as.matrix(F4c_data)
lgd1 = HeatmapAnnotation(foo=anno_block(gp = gpar(fill = 2:4),
                                        labels = c("Ctrl","10 μg/kg/d"), 
                                        labels_gp = gpar(col = "black", fontsize = 18,fontfamily='serif',fontface='bold')))
split <- factor(c("Control", "Control","Control","High dose","High dose","High dose"), levels=c("Control","High dose"))
ht=Heatmap(m.F4c_data,cluster_columns = FALSE,cluster_rows = TRUE,
           col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
           heatmap_legend_param =list(at=c(-1.5,0,1.5),title="Relative abundance",direction = "horizontal",
                                      labels_gp=gpar(fontfamily = "serif",fontsize=16),title_gp=gpar(fontface='bold',fontfamily = "serif",fontsize=18)),
           top_annotation = lgd1,show_column_names = FALSE,
           column_split =split,
           cluster_column_slices = FALSE,column_title = NULL,row_names_gp = gpar(fontfamily = "serif",fontsize=18))
ht
#png("Figure4C.jpg", width=5, height=8, units="in", res=300)
pdf('figure4c.pdf',width=5,height=6)
draw(ht,  heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()

# Figure6C
down.genes<-read.csv('./FHY/DownGenes.csv',header = TRUE)
keyvals <- ifelse(
  down.genes$Log2FC < -1, 'blue',
  ifelse(down.genes$Log2FC > 1, 'red',
         'gray39'))
keyvals[is.na(keyvals)] <- 'gray39'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'gray49'] <- 'mid'
names(keyvals)[keyvals == 'blue'] <- 'low'
keyvals.shape <- ifelse(
  down.genes$Compound =='PFO4DA', 0,
  ifelse(down.genes$Compound =='HFPO-DA', 1,
         ifelse(down.genes$Compound=='HFPO-TA',2,3)))
keyvals.shape[is.na(keyvals.shape)] <- 3
names(keyvals.shape)[keyvals.shape == 3] <- 'PFOA'
names(keyvals.shape)[keyvals.shape == 2] <- 'HFPO-TA'
names(keyvals.shape)[keyvals.shape == 1] <- 'HFPO-DA'
names(keyvals.shape)[keyvals.shape == 0] <- 'PFO4DA'
EnhancedVolcano(down.genes,
                lab =down.genes$GeneName, 
                x = 'Log2FC',
                y = 'Pvalue',
                title = '',
                selectLab = NULL,
                axisLabSize = 10,
                titleLabSize = 12,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 3,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'left',
                subtitle =NULL,
                caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                captionLabSize = 9,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = TRUE,
                shapeCustom = keyvals.shape
)+
  theme(text = element_text(family = "serif"))+
  ggplot2::coord_cartesian(xlim=c(-5,5),ylim = c(0,8))
ggsave("Figure6C.png", width = 6, height =3.7, units = c("in"), dpi = 300)
# Figure 6C.1
down.genes.1<-down.genes %>%
  filter(Compound=='PFO4DA')
keyvals <- ifelse(
  down.genes.1$Log2FC < -1, 'blue',
  ifelse(down.genes.1$Log2FC > 1, 'red',
         'gray39'))
keyvals[is.na(keyvals)] <- 'gray39'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'gray49'] <- 'mid'
names(keyvals)[keyvals == 'blue'] <- 'low'
pdf('figure6c1_1.pdf',width=4,height=3.3)
EnhancedVolcano(down.genes.1,
                lab =down.genes.1$GeneName, 
                x = 'Log2FC',
                y = 'Pvalue',
                title = 'PFO4DA',
                selectLab = c('Lpin1','Ppp1r3b','Tcim','Fkbp5',
                              'Mthfr','Tsc22d3','Tat','Slc45a3',
                              'Arrdc2','Slco1a4','Igfbp1','Zbtb16',
                              'Srebf1','Pctp','Krt23','Osbpl3','Cd36','Clstn3β'),
                axisLabSize = 12,
                titleLabSize = 14,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 0.9,
                labSize = 4,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'none',
                subtitle =NULL,
                caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                captionLabSize = 11,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = TRUE
)+
  theme(text = element_text(family = "serif"))+
  ggplot2::coord_cartesian(ylim=c(0,7.5))
dev.off()
#ggsave("Figure6C1.png", width = 4, height =3.3, units = c("in"), dpi = 300)
# Figure 6C.2
down.genes.2<-down.genes %>%
  filter(Compound=='HFPO-DA')
keyvals <- ifelse(
  down.genes.2$Log2FC < -1, 'blue',
  ifelse(down.genes.2$Log2FC > 1, 'red',
         'gray39'))
keyvals[is.na(keyvals)] <- 'gray39'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'gray49'] <- 'mid'
names(keyvals)[keyvals == 'blue'] <- 'low'
pdf('figure6c2_1.pdf',width=4,height=3.3)
EnhancedVolcano(down.genes.2,
                lab =down.genes.2$GeneName, 
                x = 'Log2FC',
                y = 'Pvalue',
                title = 'HFPO-DA',
                selectLab = c('Lpin1','Ppp1r3b','Tcim','Fkbp5',
                              'Mthfr','Tsc22d3','Tat','Slc45a3',
                              'Arrdc2','Slco1a4','Igfbp1','Zbtb16',
                              'Srebf1','Pctp','Krt23','Osbpl3','Cd36','Clstn3β'),
                axisLabSize = 12,
                titleLabSize = 14,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 0.9,
                labSize = 4,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'none',
                subtitle =NULL,
                caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                captionLabSize = 11,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = TRUE
)+
  theme(text = element_text(family = "serif"))+
  ggplot2::coord_cartesian(ylim=c(0,7.5))
dev.off()
#ggsave("Figure6C2.png", width = 4, height =3.3, units = c("in"), dpi = 300)
# Figure 6C.3
down.genes.3<-down.genes %>%
  filter(Compound=='HFPO-TA')
keyvals <- ifelse(
  down.genes.3$Log2FC < -1, 'blue',
  ifelse(down.genes.3$Log2FC > 1, 'red',
         'gray39'))
keyvals[is.na(keyvals)] <- 'gray39'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'gray49'] <- 'mid'
names(keyvals)[keyvals == 'blue'] <- 'low'
pdf('figure6c3_1.pdf',width=4,height=3.3)
EnhancedVolcano(down.genes.3,
                lab =down.genes.3$GeneName, 
                x = 'Log2FC',
                y = 'Pvalue',
                title = 'HFPO-TA',
                selectLab = c('Lpin1','Ppp1r3b','Tcim','Fkbp5',
                              'Mthfr','Tsc22d3','Tat','Slc45a3',
                              'Arrdc2','Slco1a4','Igfbp1','Zbtb16',
                              'Srebf1','Pctp','Krt23','Osbpl3','Cd36','Clstn3β'),
                axisLabSize = 12,
                titleLabSize = 14,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 0.9,
                labSize = 4,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'none',
                subtitle =NULL,
                caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                captionLabSize = 11,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = TRUE
)+
  theme(text = element_text(family = "serif"))+
  ggplot2::coord_cartesian(ylim=c(0,7.5))
dev.off()
#ggsave("Figure6C3.png", width = 4, height =3.3, units = c("in"), dpi = 300)
# Figure 6C.4
down.genes.4<-down.genes %>%
  filter(Compound=='PFOA')
keyvals <- ifelse(
  down.genes.4$Log2FC < -1, 'blue',
  ifelse(down.genes.4$Log2FC > 1, 'red',
         'gray39'))
keyvals[is.na(keyvals)] <- 'gray39'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'gray49'] <- 'mid'
names(keyvals)[keyvals == 'blue'] <- 'low'
pdf('figure6c4_1.pdf',width=4,height=3.3)
EnhancedVolcano(down.genes.4,
                lab =down.genes.4$GeneName, 
                x = 'Log2FC',
                y = 'Pvalue',
                title = 'PFOA',
                selectLab = c('Lpin1','Ppp1r3b','Tcim','Fkbp5',
                              'Mthfr','Tsc22d3','Tat','Slc45a3',
                              'Arrdc2','Slco1a4','Igfbp1','Zbtb16',
                              'Srebf1','Pctp','Krt23','Osbpl3','Cd36','Clstn3β'),
                axisLabSize = 12,
                titleLabSize = 14,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 0.9,
                labSize = 4,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'none',
                subtitle =NULL,
                caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                captionLabSize = 11,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = TRUE
)+
  theme(text = element_text(family = "serif"))+
  ggplot2::coord_cartesian(ylim=c(0,7.5))
dev.off()
#ggsave("Figure6C4.png", width = 4, height =3.3, units = c("in"), dpi = 300)
