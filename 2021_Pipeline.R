##################################################
## Schober FA & Atanassov I et al.:             ##
## A quantitative map of mitochondrial protein  ##
## phosphorylation upon loss of OXPHOS subunits ##
##################################################

################
##-- Header-- ##
################

library(ggplot2)
library(ggrepel)
library(KEGGREST)
library(pathview)
library(UniProt.ws)
library(ggridges)
library(forcats)
library(ggpubr)
library(ggseqlogo)
library(stringr)
library(limma)
library(tidyverse)
library(Hmisc)
library(cowplot)
library(org.Dm.eg.db)
library(ggalluvial)
library(pheatmap)
library(gplots)
library(gridExtra)

# Additional functions and libraries
load("Libraries/UniProt.ws.Drosophila")
source("Libraries/2021_GSEA-functions.R")
source("Libraries/2021_Plotting-functions.R")

##################
##-- Figure 1 --##
##################

# Read basic data table and meta information
info <- read.table("Data/Figure-1/info.txt", sep = "\t", header = T)
total <- read.table("Data/Figure-1/proteinGroups.txt", sep = "\t", header = T) %>%
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+") %>%
  filter(Only.identified.by.site != "+")

## -- Figure 1A -- ##
# Get the intensity columns for SILAC values
colnames <- sub("Intensity.*\\.", "", colnames(total)) 

# Which files correspond to this experiment? Narrow down dataset
files <- info[info$Experiment=="QC.TIMESERIES",]$Tube 
ts_larvae <- total[,colnames %in% files]
rownames(ts_larvae) <- total$Protein.IDs

# Data wrangling
ts_larvae %>%
  rownames_to_column(var = "ID") %>%
  gather(fraction, intensity, -ID) %>%
  mutate(sample = as.numeric(
    as.character(
      str_replace(string = str_replace(string = fraction, pattern = "L|H\\.", replacement = ""),
                  pattern = "Intensity\\.",
                  replacement = "")))) %>%
  mutate(label = str_replace_all(string = fraction, pattern = "Intensity\\.|\\.|\\d+", replacement = "")) %>%
  mutate(label = ifelse(label == "", "total", label)) %>%
  filter(label != "L") %>%
  dplyr::select(-fraction) %>%
  spread(label, intensity) %>%
  drop_na(total) %>%
  filter(total > 0) %>%
  mutate(fraction = H/total) %>%
  left_join(info, by = c("sample" = "Tube")) -> ts_larvae_plot

cairo_pdf("Figures/Figure-1A.pdf", 4, 3)
ggplot(data = ts_larvae_plot, aes(x = Condition, y = 100*fraction))+
  stat_summary(aes(shape = Replicate),
               fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom = "pointrange",  size = 0.6,
               position = position_dodge(0.8))+
  theme_classic()+
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0,108))+
  geom_hline(yintercept = 95, lty = "dotted")
dev.off()

ts_larvae_plot %>%
  group_by(Condition) %>%
  summarise(mean = mean(fraction), sd = sd(fraction))

## -- Figure 1B -- ##
# Get the intensity columns for SILAC values
colnames <- sub("Intensity.*\\.", "", colnames(total)) 

# Which files correspond to this experiment?
files <- info[info$Experiment=="QC.TIMESERIES.FLIES",]$Tube 
ts_f <- total[,colnames %in% files]
rownames(ts_f) <- total$Protein.IDs

# Data wrangling
ts_f %>%
  rownames_to_column(var = "ID") %>%
  gather(fraction, intensity, -ID) %>%
  mutate(sample = as.numeric(
    as.character(
      str_replace(string = str_replace(string = fraction, pattern = "L|H\\.", replacement = ""),
                  pattern = "Intensity\\.",
                  replacement = "")))) %>%
  mutate(label = str_replace_all(string = fraction, pattern = "Intensity\\.|\\.|\\d+", replacement = "")) %>%
  mutate(label = ifelse(label == "", "total", label)) %>%
  filter(label != "L") %>%
  dplyr::select(-fraction) %>%
  spread(label, intensity) %>%
  drop_na(total) %>%
  filter(total > 0) %>%
  mutate(fraction = H/total) %>%
  left_join(info, by = c("sample" = "Tube")) %>%
  mutate(SYMBOL = mapIds(org.Dm.eg.db, 
                          keys=as.character(sub(";.*", "", ID)), 
                          column="SYMBOL", 
                          keytype="UNIPROT",
                          multiVals="first"))-> ts_f_plot

ts_f_plot %>%
  group_by(Condition) %>%
  summarise(mean = mean(fraction), sd = sd(fraction)) -> ts_f_plot_quant

# Plotting
pdf("Figures/Figure-1B.pdf")
ggplot(ts_f_plot, aes(x = fraction, y = fct_rev(Condition))) + 
  geom_density_ridges_gradient(aes(fill = ..x..), scale = 3, size = 0.5, panel_scaling = T)+
  scale_fill_gradientn(colours = c("white", "#D54036"), name = "Labelling efficiency")+
  theme_ridges(font_size = 13, grid = T) + theme(axis.title.y = element_blank())+
  xlab("Labelling efficiency")+
  scale_x_continuous(breaks=seq(from = 0, to = 1, by = 0.2))+
  ggtitle("Figure 1B: Timecourse labelling efficiency in f")
dev.off()

## -- Figure 1C -- ##
# Which categories are overrepresented in the four bins [0-0.25], [0.25-0.5], [0.5-0.75], [0.75-1.0]?

ts_f_plot %>%
  group_by(Condition, ID) %>%
  summarise(mean = mean(fraction)) %>%
  filter(Condition == "d14") %>%
  mutate(ID_short = sub(";.*", "", ID)) -> ts_f_ora_d14

entrez.tmp <- UniProt.ws::select(up, ts_f_ora_d14$ID_short, "ENTREZ_GENE")
left_join(ts_f_ora_d14, entrez.tmp, by = c("ID_short" = "UNIPROTKB")) -> ts_f_ora

# Define background
write.table(ts_f_ora$ENTREZ_GENE,
            quote = F, row.names = F, col.names = F, file = "Data/Figure-1/GSEA/ts_f_total.txt")

# Write entrez ID lists for input to Webgestalt
write.table(ts_f_ora[ts_f_ora$mean >= 0 & ts_f_ora$mean < 0.25,]$ENTREZ_GENE,
            quote = F, row.names = F, col.names = F, file = "Data/Figure-1/GSEA/ts_f_0-0.25.txt")
write.table(ts_f_ora[ts_f_ora$mean >= 0.25 & ts_f_ora$mean < 0.5,]$ENTREZ_GENE,
            quote = F, row.names = F, col.names = F, file = "Data/Figure-1/GSEA/ts_f_0.25-0.5.txt")
write.table(ts_f_ora[ts_f_ora$mean >= 0.5 & ts_f_ora$mean < 0.75,]$ENTREZ_GENE,
            quote = F, row.names = F, col.names = F, file = "Data/Figure-1/GSEA/ts_f_0.5-0.75.txt")
write.table(ts_f_ora[ts_f_ora$mean >= 0.75 & ts_f_ora$mean < 1,]$ENTREZ_GENE,
            quote = F, row.names = F, col.names = F, file = "Data/Figure-1/GSEA/ts_f_0.75-1.txt")

# Read in categories and annotate the proteins detected in d14 flies
ts.f_kegg_q1 <- read.table("Data/Figure-1/GSEA/Timeseries_flies_webgestalt_KEGG_0-0.25/enrichment_results_wg_result1610357770.txt", header = T, sep = "\t") %>%
  filter(FDR < 0.05)
ts.f_kegg_q2 <- read.table("Data/Figure-1/GSEA/Timeseries_flies_webgestalt_KEGG_0.25-0.5/enrichment_results_wg_result1610357778.txt", header = T, sep = "\t") %>%
  filter(FDR < 0.05)
ts.f_kegg_q3<- read.table("Data/Figure-1/GSEA/Timeseries_flies_webgestalt_KEGG_0.5-0.75/enrichment_results_wg_result1610357785.txt", header = T, sep = "\t") %>%
  filter(FDR < 0.05)
ts.f_kegg_q4 <- read.table("Data/Figure-1/GSEA/Timeseries_flies_webgestalt_KEGG_0.75-1/enrichment_results_wg_result1610357790.txt", header = T, sep = "\t") %>%
  filter(FDR < 0.05)

# Data wrangle
ts.f_kegg <- rbind(ts.f_kegg_q1, ts.f_kegg_q2, ts.f_kegg_q3, ts.f_kegg_q4)

ts_f_plot %>%
  group_by(Condition, ID) %>%
  summarise(mean = mean(fraction)) %>%
  mutate(ID_short = sub(";.*", "", ID)) %>%
  left_join(entrez.tmp, by = c("ID_short" = "UNIPROTKB")) -> ts_f_ora_all

entrez_gene_kegg_annot = "Libraries/dmelanogaster_pathway_KEGG_entrezgene.gmt"
entrez_gene_kegg_annot <- read_lines(entrez_gene_kegg_annot)

entrez_gene_kegg_annot %>%
  map(str_split,pattern = "\\\t") %>% 
  map(unlist) %>% 
  map(1) -> kegg_ids_list

entrez_gene_kegg_annot %>%
  map2(.x = .,
       .y = kegg_ids_list,
       .f = kegg_ids_tab_to_long) %>% 
  map(as_tibble) %>% 
  map(separate_rows, value, sep = "\\\t") %>% 
  map2(.x = .,
       .y = kegg_ids_list,
       .f =  add_kegg_id) %>% 
  bind_rows() %>% 
  dplyr::rename(geneSet = kegg_id,
                Id = value) %>% 
  dplyr::select(geneSet, Id) -> entrez_gene_kegg_annot_long

ts_f_ora_all %>%
  mutate(ENTREZ_GENE = str_replace_all(ENTREZ_GENE, " ", "")) %>%
  left_join(entrez_gene_kegg_annot_long, by = c("ENTREZ_GENE" = "Id")) %>%
  left_join(ts.f_kegg) %>%
  drop_na(description) -> ts_f_ora_plot

ts_f_ora_plot %>%
  group_by(Condition, description) %>%
  summarise(median = median(mean)) -> ts_f_ora_stats

# Define colors
colfunc <- colorRampPalette(c("white", "#D54036"))
color <- data.frame(median = round(seq(0,1, by = 0.01),2), color = colfunc(101))

ts_f_ora_plot %>%
  group_by(description) %>%
  summarise(median = median(mean)) %>%
  arrange(median) -> col_sum

ts_f_ora_plot %>%
  mutate(description = fct_relevel(description, col_sum$description)) %>%
  group_by(Condition, description) %>%
  summarise(median = median(mean)) %>%
  mutate(median = round(median, 2)) %>%
  left_join(color) %>%
  arrange(description, Condition) -> col_sum_cond

cairo_pdf("Figures/Figure-1C.pdf", 15, 4)
ts_f_ora_plot %>%
  mutate(description = fct_relevel(description, col_sum$description)) %>%
  ggplot(aes (x = Condition, y = mean))+
  geom_boxplot(fill = col_sum_cond$color)+
  theme_classic()+
  scale_y_continuous(limits = c(0,1))+
  geom_hline(yintercept = 0.25, lty = "dotted") + 
  geom_hline(yintercept = 0.5, lty = "dotted") + 
  geom_hline(yintercept = 0.75, lty = "dotted") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  facet_wrap(.~description, nrow = 1)
dev.off()

ts_f_plot %>%
  mutate(ENTREZ_GENE = mapIds(org.Dm.eg.db, 
                         keys=as.character(sub(";.*", "", ID)), 
                         column="ENTREZID", 
                         keytype="UNIPROT",
                         multiVals="first")) %>%
  mutate(ENTREZ_GENE = str_replace_all(ENTREZ_GENE, " ", "")) %>%
  left_join(entrez_gene_kegg_annot_long, by = c("ENTREZ_GENE" = "Id")) %>%
  left_join(ts.f_kegg) -> ts_f_table_s2
  

write.table(ts_f_table_s2,file = "SuppTables/Fig1C.txt", row.names = F, quote = F, sep = "\t")

###################
##-- Figure S1 --##
###################

## -- Figure S1A and B -- ##
## -- Show conversion into Proline -- ##
mod.spec.pept <- read_tsv("Data/Figure-S1/Proline/modificationSpecificPeptides.txt")
mod.spec.pept$containsK <- grepl("K", mod.spec.pept$Sequence)
mod.spec.pept$containsKR <- grepl("R|K", mod.spec.pept$Sequence)
mod.spec.pept$containsP <- grepl("P", mod.spec.pept$Sequence)

mod.spec.pept %>%
  filter(containsKR == TRUE) %>%
  mutate(`Modification state` = ifelse(Modifications %in% "Unmodified","Unmodified",
                                       ifelse(grepl("Pro6",Modifications),"Pro6 Modified","Arg / Lys Modified"))) %>% 
  group_by(`Raw file`,`Modification state`) %>% 
  summarise(Count = n()) %>% 
  ungroup() %>% 
  group_by(`Raw file`) %>% 
  mutate(`Normalized count [%]` = 100*(`Count`/sum(`Count`))) %>% 
  filter(grepl("31|37",`Raw file`)) %>% 
  ggplot(aes(x = `Raw file`, y = `Normalized count [%]`, fill = `Modification state`))+
  geom_bar(stat = "identity", color = "black")+
  geom_text(aes(label = round(`Normalized count [%]`,2)),position = position_stack(vjust = 0.5))+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
ggsave(filename = "Figures/Figure-S1A-A.pdf")

mod.spec.pept %>%
  filter(containsK == TRUE) %>%
  mutate(`Modification state` = ifelse(Modifications %in% "Unmodified","Unmodified",
                                       ifelse(grepl("Pro6",Modifications),"Pro6 Modified","Arg / Lys Modified"))) %>% 
  group_by(`Raw file`,`Modification state`) %>% 
  summarise(Count = n()) %>% 
  ungroup() %>% 
  group_by(`Raw file`) %>% 
  mutate(`Normalized count [%]` = 100*(`Count`/sum(`Count`))) %>% 
  filter(grepl("101|108",`Raw file`)) %>% 
  ggplot(aes(x = `Raw file`, y = `Normalized count [%]`, fill = `Modification state`))+
  geom_bar(stat = "identity", color = "black")+
  geom_text(aes(label = round(`Normalized count [%]`,2)),position = position_stack(vjust = 0.5))+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
ggsave(filename = "Figures/Figure-S1A-B.pdf")

mod.spec.pept %>%
  filter(containsP == TRUE) %>%
  filter(containsKR == TRUE) %>%
  mutate(`Modification state` = ifelse(Modifications %in% "Unmodified","Unmodified",
                                       ifelse(grepl("Pro6",Modifications),"Pro6 Modified","Arg / Lys Modified"))) %>% 
  group_by(`Raw file`,`Modification state`) %>% 
  summarise(Count = n()) %>% 
  ungroup() %>% 
  group_by(`Raw file`) %>% 
  mutate(`Normalized count [%]` = 100*(`Count`/sum(`Count`))) %>% 
  filter(grepl("31|37",`Raw file`)) %>% 
  ggplot(aes(x = `Raw file`, y = `Normalized count [%]`, fill = `Modification state`))+
  geom_bar(stat = "identity", color = "black")+
  geom_text(aes(label = round(`Normalized count [%]`,2)),position = position_stack(vjust = 0.5))+
  scale_fill_brewer(palette = "Paired")+
  theme_classic()+
ggsave(filename = "Figures/Figure-S1B_A.pdf")

mod.spec.pept %>%
  filter(containsP == TRUE) %>%
  filter(containsK == TRUE) %>%
  mutate(`Modification state` = ifelse(Modifications %in% "Unmodified","Unmodified",
                                       ifelse(grepl("Pro6",Modifications),"Pro6 Modified","Arg / Lys Modified"))) %>% 
  group_by(`Raw file`,`Modification state`) %>% 
  summarise(Count = n()) %>% 
  ungroup() %>% 
  group_by(`Raw file`) %>% 
  mutate(`Normalized count [%]` = 100*(`Count`/sum(`Count`))) %>% 
  filter(grepl("101|108",`Raw file`)) %>% 
  ggplot(aes(x = `Raw file`, y = `Normalized count [%]`, fill = `Modification state`))+
  geom_bar(stat = "identity", color = "black")+
  geom_text(aes(label = round(`Normalized count [%]`,2)),position = position_stack(vjust = 0.5))+
  scale_fill_brewer(palette = "Paired")+
  theme_classic()+
ggsave(filename = "Figures/Figure-S1B-B.pdf")
  
## -- Figure S1C -- ##
files <- info[grepl("SILAF.PRECISION|DOUBLE.PRECISION",info$Experiment),]$Tube
colnames <- sub("Ratio.H.L.", "", colnames(total)[grepl("Ratio.H.L.\\d", colnames(total))])
total[,grepl("Ratio.H.L.\\d", colnames(total))] -> KR_ratio_ds
KR_ratio_ds <- KR_ratio_ds[,colnames %in% files]
rownames(KR_ratio_ds) <- total$Protein.IDs

cairo_pdf("Figures/Figure-S1C.pdf", 4.5, 3)
KR_ratio_ds %>%
  rownames_to_column("ID") %>%
  gather(Sample, Ratio, 2:8) %>%
  drop_na(Ratio) %>%
  mutate(Sample = as.numeric(as.character(str_replace_all(Sample, "Ratio.H.L.", "")))) %>%
  left_join(info %>% filter(!grepl("LFQ", info$Experiment)), by = c("Sample" = "Tube")) %>%
  ggplot(aes(x = Experiment, y = log2(Ratio), fill = Experiment))+
    geom_violin()+
  stat_summary(fun=median, geom="point", size=2, color="black", shape = 15)+
  theme_classic()+
  geom_hline(yintercept = 0, lty = "dotted")+
  scale_y_continuous(breaks = seq(-4,3, by = 1))
dev.off()

## -- Figure S1D -- ##
files <- info[grepl("SILAF.PRECISION|DOUBLE.PRECISION",info$Experiment),]$Tube
colnames <- sub("Intensity.L.", "", colnames(total)[grepl("Intensity.L.", colnames(total))])
total[,grepl("Intensity.L.", colnames(total))] -> KR_int_ds
KR_int_ds <- KR_int_ds[,colnames %in% files]
rownames(KR_int_ds) <- total$Protein.IDs
  
cairo_pdf("Figures/Figure-S1D.pdf", 4.5, 3)
KR_int_ds %>%
  rownames_to_column("ID") %>%
  gather(Sample, Intensity, 2:8) %>%
  drop_na(Intensity) %>%
  filter(Intensity > 0) %>%
  mutate(Sample = as.numeric(as.character(str_replace_all(Sample, "Intensity.L.", "")))) %>%
  left_join(info %>% filter(!grepl("LFQ", info$Experiment)), by = c("Sample" = "Tube")) -> KR_int_ds_plot

KR_int_ds_plot %>%
    ggplot(aes(x = Experiment, y = log2(Intensity), fill = Experiment))+
    geom_violin()+
    stat_summary(fun=median, geom="point", size=2, color="black", shape = 15)+
    theme_classic()+
    geom_hline(yintercept = 0, lty = "dotted")
dev.off()

KR_int_ds_plot %>%
  group_by(Experiment) %>%
  summarise(mean = mean(Intensity))

## -- Figure S1E -- ##
## LFQ standard deviations
files.current <- info[info$Experiment=="QC.LFQ.SILAFL",][,"Tube"]
colnames <- sub("LFQ\\.intensity\\.L\\.", "", colnames(total))
LFQ <- total[,colnames %in% files.current]
LFQ[LFQ == "NaN" | LFQ == 0] <- NA

### Normalize to 1
LFQ.mean <- apply(LFQ, 1, function(x) mean(x, na.rm = T))
LFQ <- LFQ/LFQ.mean

LFQ <- LFQ[!rowSums(is.na(LFQ)) > 1,]
LFQ.sd <- apply(LFQ, 1, function(x) sd(x, na.rm = T))

## SILAF standard deviations
files.current <- info[info$Experiment=="QC.SILAF.PRECISION",][,"Tube"]
colnames <- sub("Ratio\\.H\\.L\\.normalized\\.", "", colnames(total))
SILAF <- total[,colnames %in% files.current]
SILAF[SILAF == 0 | SILAF== "NaN"] <- NA
SILAF <- SILAF[!rowSums(is.na(SILAF)) > 1,]

SILAF.sd <- apply(SILAF, 1, function(x) sd(x, na.rm = T))

## Plot density
precision <- data.frame(SD = c(LFQ.sd, SILAF.sd),
                        dataset = c(rep("LFQ", length(LFQ.sd)), rep("SILAF", length(SILAF.sd))))

cairo_pdf("Figures/Figure-S1E.pdf")
ggplot(precision, aes(x = SD, fill = dataset)) +
  geom_density(alpha = 0.5)+
  coord_cartesian(xlim = c(0,1.0))+
  theme_classic()+
  scale_fill_manual(values = c("grey", "firebrick3"))
dev.off()

## -- Figure S1F -- ##
### (1) Read in data >>>
files <- info[grepl("QC.SILAF.PRECISION",info$Experiment),]$Tube
colnames <- sub("Intensity.|Intensity.H.|Intensity.L.", "", colnames(total)[grepl("Intensity.", colnames(total))])
total[,grepl("Intensity.", colnames(total))] -> precision.SILAF
precision.SILAF <- log2(precision.SILAF[,colnames %in% files])
precision.SILAF[abs(precision.SILAF) == Inf] <- NA

### (2) Split the samples into elements of a list. Necessary to generate seperate scatterplots >>>

file.counter <- length(unique(sub(".*\\.","",colnames(precision.SILAF)))) # How many samples?

## Iterate column groups for one sample into one element of a list each. There are file.counter = i samples >>
for (i in c(1:file.counter)){
  if(i==1){
    precision.SILAF.samples <- list() # Empty list in first round
  }
  data.current <- precision.SILAF[,(i*3-2):(i*3)] # Narrow down to the current dataset
  data.current <- data.current[data.current[,1]!=0,] # Exclude all rows with a total count of 0
  sample.current <- sub(".*\\.","",colnames(data.current)[1]) # What's the samplename?
  precision.SILAF.samples[[sample.current]] <- data.current # Give the listelement the name of the current sample
  rm(data.current, sample.current)
}

### (3) Generates one plot per element in the list >>>
cairo_pdf("Figures/Figure-S1F.pdf", 6,9)
par(mfrow = c(4,3))
for (i in c(1:file.counter)){
  panel.smooth(precision.SILAF.samples[[i]][,2], precision.SILAF.samples[[i]][,3])
  panel.hist(precision.SILAF.samples[[i]][,2], c(15,20,25,30,35,40))
  panel.hist(precision.SILAF.samples[[i]][,3], c(15,20,25,30,35,40))
}
dev.off()

## -- Figure S1G -- ##
files <- info[grepl("QC.LFQ.SILAFL",info$Experiment),]$Tube
total %>%
  dplyr::select(contains("LFQ.intensity.L")) %>%
  dplyr::select(contains(as.character(files))) -> precision.LFQ.holidic

precision.LFQ.holidic[precision.LFQ.holidic == 0] <- NA

cairo_pdf("Figures/Figure-S1G.pdf")
pairs(precision.LFQ.holidic,
      lower.panel = panel.smooth.add,
      diag.panel = panel.hist.add,
      upper.panel = NULL,
      main = "S1D: LFQ samples on holidic food")
dev.off()

###################
##-- Figure S2 --##
###################

## -- Figure S2A -- ##
ts_larvae[,!grepl(".*[LH].*", colnames(ts_larvae))] %>%
  gather("sample", "Intensity") %>%
  filter(Intensity > 0) %>%
  mutate(sample = as.numeric(str_replace_all(sample, "Intensity\\.", ""))) %>%
  filter(sample > 101) %>%
  arrange(sample, Intensity) -> dynamic

rank_table <- as.vector(table(dynamic$sample))
for(i in 1:length(rank_table)){
  if(i == 1){rank_vector <- c()}
  rank_vector <- c(rank_vector, (1:rank_table[i]))
}

dynamic %>%
  mutate(rank = rank_vector) -> dynamic_plot

cairo_pdf("Figures/Figure-S2A.pdf", 5, 5)
ggplot(data = dynamic_plot, aes (x = rank, y = log2(Intensity), col = as.character(sample)))+
  geom_point()+
  scale_color_brewer(palette = "Reds")+
  theme_bw()
dev.off()

## -- Figure S2B -- ##
# Read data
read.table("Data/Figure-S2/proteinGroups.txt", sep = "\t", header = T) %>%
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+") %>%
  filter(Only.identified.by.site != "+") -> flysex

flysex[flysex == "NaN"] <- NA

# Entrez ID mapping
entrez <- UniProt.ws::select(up, sub(";.*", "", flysex$Protein.IDs), "ENTREZ_GENE")
entrez <- entrez[!duplicated(entrez$UNIPROTKB),]
flysex$entrez <- entrez$ENTREZ_GENE[match(entrez$UNIPROTKB, sub(";.*", "", flysex$Protein.IDs))]

flysex %>%
  dplyr::select(matches("normalized.SILAF|Gene.names|ENTREZ|Protein")) %>%
  dplyr::rename(Rep1 = Ratio.H.L.normalized.SILAF.F.vs.M.Fod,
                Rep2 = Ratio.H.L.normalized.SILAF.F.vs.M.Rev) %>%
  mutate(Rep1 = log2(Rep1), Rep2 = -log2(Rep2)) -> flysex_plot

rowsums <- rowSums(is.na(flysex_plot %>% dplyr::select(contains("Rep"))))
proteins.n <- as.character(as.numeric((dim(flysex_plot %>% filter(rowsums < 2))[1])))

print(paste("Depth of female/male comparison:", dim(flysex_plot %>% filter(rowsums < 2))[1]))

flysex_plot %>%
  filter(rowSums(is.na(flysex_plot %>% dplyr::select(contains("Rep")))) < 2) %>%
  rowwise %>%
  mutate(mean = mean(c(Rep1, Rep2), na.rm = T)) %>%
  write.table(., file = "SuppTables/FigS2.txt", row.names = F, sep = "\t", quote = F)

#Calculate R:
round(cor(flysex_plot$Rep1, flysex_plot$Rep2, use = "complete.obs"), digits=4)

cairo_pdf("Figures/Figure-S2B.pdf", 6, 6)
ggplot()+
  geom_point(data = flysex_plot %>% filter(abs(Rep1) <= 4 | abs(Rep2) <= 4), aes(x = Rep1, y = Rep2), color="black", alpha = 0.3) +
  geom_point(data = flysex_plot %>% filter(abs(Rep1) > 4 & abs(Rep2) > 4), aes(x = Rep1, y = Rep2), color = "firebrick3")+
  geom_abline(intercept = 0, slope = 1, lty = "dotted")+
  geom_hline(yintercept = 0, col = "black")+
  geom_vline(xintercept = 0, col = "black")+
  geom_point(alpha = 0.2)+
  coord_cartesian(xlim = c(-8,8), ylim = c(-8,8))+
  geom_text_repel(
    data = flysex_plot %>% filter(abs(Rep1) > 4 & abs(Rep2) > 4),
    aes(x = Rep1, y = Rep2, label = Gene.names),
    size = 1)+
  theme_classic() 
dev.off()

## -- Figure S2C -- ##
## Generate the table to upload to WebGestalt
flysex_gsea_in <- data.frame(ENTREZ_GENE = flysex_plot$entrez, mean_sthlm = apply(flysex_plot[,6:7], 1, function(x) mean(x, na.rm = T)))
write.table(flysex_gsea_in, "Data/Figure-S2/GSEA/GSEA.Flysex.deep.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

## GSEA plots: Barchart and volcano plot
read.table("Data/Figure-S2/GSEA/Project_wg_result1610481328/enrichment_results_wg_result1610481328.txt", sep = "\t", header = T) %>%
  mutate(description = fct_reorder(description, normalizedEnrichmentScore)) -> flysex_gsea_plot

cairo_pdf("Figures/Figure-S2C.pdf", 6, 7)
ggplot(flysex_gsea_plot, aes(x=normalizedEnrichmentScore, y=description, label=description)) + 
  geom_point(stat="identity", size=4)  +
  geom_segment(aes(y = description, 
                   x = normalizedEnrichmentScore, 
                   yend = description, 
                   xend = 0))+
  theme_classic()+
  geom_vline(xintercept = 0, lty = "dotted")
dev.off()

## -- Figure S2D -- ##
# Import the Sury 2010 et al. dataset (doi: 10.1074/mcp.M110.000323)
read.table("Data/Figure-S2/IA_FS_our_SILAF_published_SILAF/proteinGroups.txt", header = T, sep = "\t") %>%
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+") %>%
  filter(Only.identified.by.site != "+") %>%
  dplyr::rename(w1118.F = Ratio.H.L.normalized.w1118_female_standard, w1118.M = Ratio.H.L.normalized.w1118_male_standard) %>%
  dplyr::select(Protein.IDs, Gene.names, w1118.F, w1118.M) %>%
  mutate(Ratio_old = log2(w1118.M/w1118.F)) %>%
  drop_na(Ratio_old) %>%
  mutate(Protein.IDs_simple = sub(";.*", "", Protein.IDs)) -> flysex_old

entrez <- UniProt.ws::select(up, flysex_old$Protein.IDs_simple, "ENTREZ_GENE") #Maps entrez IDs to uniprot IDs
entrez <- entrez[!duplicated(entrez$UNIPROTKB),]

flysex_old <- left_join(flysex_old, entrez, by = c("Protein.IDs_simple" = "UNIPROTKB")) %>% drop_na(ENTREZ_GENE)

flysex_comparative <- left_join(flysex_old, flysex_gsea_in %>% drop_na(mean_sthlm))

#Calculate R:
  round(cor(flysex_comparative$Ratio_old, flysex_comparative$mean_sthlm, use = "complete.obs"), digits=4)

  ## -- Figure S2D -- ##
cairo_pdf("Figures/Figure-S2D.pdf")
ggplot()+
  geom_point(data = flysex_comparative %>% filter(abs(Ratio_old) > 2 | abs(mean_sthlm) > 2),
             aes (x = mean_sthlm, y = Ratio_old), color = "red") +
  geom_point(data = flysex_comparative %>% filter(abs(Ratio_old) <= 2 & abs(mean_sthlm) <= 2),
             aes (x = mean_sthlm, y = Ratio_old), color = "black", alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1)+
  geom_text_repel(data = flysex_comparative %>% filter(abs(Ratio_old) > 2 | abs(mean_sthlm) > 2),
                  aes (x = mean_sthlm, y = Ratio_old, label = Gene.names),
                  size = 3)+
  theme_classic()+
  geom_smooth(method='lm',formula=y~x)+
  coord_cartesian(xlim = c(-10,10), ylim = c(-10,10))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)
dev.off()

########################
##-- Figure 2 and S3--##
########################

## Data wrangle
read.table("Data/Figure-2/proteinGroups.txt", sep = "\t", header = T) %>%
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+") %>%
  filter(Only.identified.by.site != "+") %>%
  dplyr::rename(Rep1 = Ratio.H.L.normalized.1, Rep2 = Ratio.H.L.normalized.2) %>%
  dplyr::select(Protein.IDs, Protein.names, Gene.names, Rep1, Rep2, contains("Intensity.L."), contains("Intensity.H.")) %>%
  filter(complete.cases(Rep1, Rep2)) %>%
  mutate(Rep1 = log2(Rep1), Rep2 = log2(Rep2)) -> holidic_proteome

## Quantities for Figure 2B
tmp <- dim(read.table("Data/Figure-2/proteinGroups.txt", sep = "\t", header = T) %>%
      filter(Potential.contaminant != "+") %>%
      filter(Reverse != "+") %>%
      filter(Only.identified.by.site != "+"))[1]
print(paste("All detected proteins:", tmp))
tmp <- dim(read.table("Data/Figure-2/proteinGroups.txt", sep = "\t", header = T) %>%
      filter(Potential.contaminant != "+") %>%
      filter(Reverse != "+") %>%
      filter(Only.identified.by.site != "+") %>%
      dplyr::rename(Rep1 = Ratio.H.L.normalized.1, Rep2 = Ratio.H.L.normalized.2) %>%
      dplyr::select(Protein.IDs, Protein.names, Gene.names, Rep1, Rep2) %>%
      mutate(NANA = paste(Rep1, Rep2)) %>%
      filter(NANA != "NaN NaN"))[1]
print(paste("Quantified at least once:", tmp))
tmp <- dim(read.table("Data/Figure-2/proteinGroups.txt", sep = "\t", header = T) %>%
             filter(Potential.contaminant != "+") %>%
             filter(Reverse != "+") %>%
             filter(Only.identified.by.site != "+") %>%
             dplyr::rename(Rep1 = Ratio.H.L.normalized.1, Rep2 = Ratio.H.L.normalized.2) %>%
             dplyr::select(Protein.IDs, Protein.names, Gene.names, Rep1, Rep2) %>%
             filter(complete.cases(Rep1, Rep2)))[1]
print(paste("Quantified at least once:", tmp))

## Limma input
hol_prot_limma <- holidic_proteome %>% column_to_rownames("Protein.IDs") %>% dplyr::select(Rep1, Rep2)

## Calculate padj with limma moderated t-test
fit <- lmFit(hol_prot_limma)
fit <- eBayes(fit)
topTable(fit, number = Inf, confint = TRUE, adjust.method = "fdr") %>%
  rownames_to_column("Protein.IDs") %>%
  dplyr::select(Protein.IDs, logFC, adj.P.Val) %>%
  right_join(holidic_proteome) -> holidic_proteome

write.table(holidic_proteome, file = "SuppTables/Fig2.txt", row.names = F, quote = F, sep = "\t")

table(holidic_proteome$adj.P.Val <= 0.01)
table(holidic_proteome$adj.P.Val <= 0.01 & holidic_proteome$logFC > 0)
table(holidic_proteome$adj.P.Val <= 0.01 & holidic_proteome$logFC < 0)

## Prepare a table for GSEA
rankedSet <- UniProt.ws::select(up, sub(";.*","", holidic_proteome$Protein.IDs), "ENTREZ_GENE")
holidic_proteome %>%
  mutate(Protein.IDs_simple = sub(";.*","", Protein.IDs)) %>%
  left_join(rankedSet, by = c("Protein.IDs_simple" = "UNIPROTKB")) %>%
  mutate(ENTREZ_GENE = as.numeric(as.character(ENTREZ_GENE))) %>%
  arrange(logFC) -> holidic_proteome

## Generate table for Webgestalt GSEA
holidic_proteome %>% 
  dplyr::select(ENTREZ_GENE, logFC) -> holidic_gsea_webgestalt
  
write.table(holidic_gsea_webgestalt, "Data/Figure-2/GSEA_proteome//Holidic_GSEA.rnk", quote = F, col.names = F, row.names = F, sep = "\t")

## Import KEGG restults and represent them
holidic_gsea_result <- read.table("Data/Figure-2/GSEA_proteome/KEGG_unlog/enrichment_results_wg_result1610619990.txt", header = T, sep = "\t") %>%
  filter(FDR < 0.05)
holidic_quant_data <- read_tsv("Data/Figure-2/GSEA_proteome/KEGG_unlog/interestingID_mappingTable_wg_result1610619990.txt")

## -- Figure 2B -- ##
cairo_pdf("Figures/Figure-2B.pdf")
ggplot() +
  geom_point(data = holidic_proteome %>% filter(adj.P.Val > 0.01), aes(x = logFC, y = -log10(adj.P.Val)),
             color = "grey30", alpha = 0.6)+
  geom_point(data = holidic_proteome %>% filter(adj.P.Val <= 0.01 & logFC < 0), aes(x = logFC, y = -log10(adj.P.Val)),
             color = "#1F80AA", alpha = 0.6)+
  geom_point(data = holidic_proteome %>% filter(adj.P.Val <= 0.01 & logFC > 0), aes(x = logFC, y = -log10(adj.P.Val)),
             color = "firebrick3", alpha = 0.6)+
  theme_classic()+
  labs(x = "log2(fold change)", y = "-log10(adjusted p-value)")+
  annotate(geom="text", x=-4, y=2.8, label= dim(holidic_proteome %>% filter(adj.P.Val <= 0.01 & logFC < 0) %>% dplyr::select(logFC))[1], color="#1F80AA")+
  annotate(geom="text", x=4, y=2.8, label= dim(holidic_proteome %>% filter(adj.P.Val <= 0.01 & logFC > 0) %>% dplyr::select(logFC))[1], color="firebrick3")+
  geom_hline(yintercept = -log10(0.01), color = "black", lty = "dotted")+
  coord_cartesian(xlim = c(-4,4))
dev.off()

## -- Figure 2C -- ##
cairo_pdf("Figures/Figure-2C.pdf", 5, 7)
plot_gsea(entrez_gene_kegg_annot = "Libraries/dmelanogaster_pathway_KEGG_entrezgene.gmt",
          quant_data = "Data/Figure-2/GSEA_proteome/KEGG/interestingID_mappingTable_wg_result1610546767.txt",
          gsea_result = "Data/Figure-2/GSEA_proteome/KEGG/enrichment_results_wg_result1610546767.txt")
dev.off()

## -- Figure S3A -- ##
holidic_phospho <- read_tsv("Data/Figure-2/Phospho (STY)Sites.txt") %>%
  filter(`Localization prob`>=0.75) %>%
  filter(!grepl("REV__", Protein)) %>%
  mutate(SYMBOL = mapIds(org.Dm.eg.db, 
                         keys=as.character(Protein), 
                         column="SYMBOL", 
                         keytype="UNIPROT",
                         multiVals="first"))
  
holidic_phospho[holidic_phospho == "NaN"] <- NA

cairo_pdf("Figures/Figure-S3A.pdf")
holidic_phospho %>%
  filter(str_length(`Sequence window`) == 31) %>%
  drop_na(`Sequence window`) %>%
ggseqlogo(data = .$`Sequence window`, method = 'bits', col_scheme='chemistry')
dev.off()

## -- Figure S3B -- ##
#Calculate R:
correl <- round(cor(log2(holidic_phospho$`Ratio H/L normalized 1`), log2(holidic_phospho$`Ratio H/L normalized 2`), use = "complete.obs"), digits=4)

cairo_pdf("Figures/Figure-S3B.pdf")
ggplot(data = holidic_phospho, aes(x = log2(`Ratio H/L normalized 1`), y = log2(`Ratio H/L normalized 2`)))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_abline(slope = 1, intercept = 0, lty = "dotted")+
  theme_classic()+
  geom_point(data = subset(holidic_phospho, abs(log2(`Ratio H/L normalized 1`)) < 2 &
                             abs(log2(`Ratio H/L normalized 2`)) < 2),
             color = "black", alpha = 0.3)+
  geom_point(data = subset(holidic_phospho, abs(log2(`Ratio H/L normalized 1`)) >= 2 &
                             abs(log2(`Ratio H/L normalized 2`)) >= 2),
             color = "firebrick3", alpha = 1)+
  geom_text_repel(
    data = subset(holidic_phospho, abs(log2(`Ratio H/L normalized 1`)) >= 2 &
                    abs(log2(`Ratio H/L normalized 2`)) >= 2),
    aes(label = SYMBOL),
    size = 5,
    box.padding = 0.5,
    point.padding = 0.5,
    max.overlaps = 15)+
  coord_fixed(ratio = 1/1)+
  annotate(geom="text", x=-4, y=3, label= paste("R = ", correl), color="black")+
  scale_x_continuous(limits = c(-5,5), breaks = seq(-5,5,by = 1), name = "Replicate 1 (log2 Holidic/Yeast food)") +
  scale_y_continuous(limits = c(-5,5), breaks = seq(-5,5,by = 1), name = "Replicate 2 (log2 Holidic/Yeast food)")
dev.off()

## -- Figure S3C -- ##
holidic_phospho[holidic_phospho == "NaN"] <- NA

phospho_occupancy <- holidic_phospho %>%
  dplyr::select(Protein, Proteins, `Positions within proteins`, `Occupancy H 1`, `Occupancy L 1`,
         `Occupancy H 2`, `Occupancy L 2`, SYMBOL, `Ratio H/L normalized 1`, `Ratio H/L normalized 2`) %>%
  mutate(Rep1_occ = `Occupancy H 1`/`Occupancy L 1`, Rep2_occ = `Occupancy H 2`/`Occupancy L 2`) %>%
  mutate(SYMBOL = paste(SYMBOL)) %>%
  rowwise() %>%
  mutate(occ_mean = log2(mean(c(Rep1_occ, Rep2_occ), na.rm =T)), ratio_mean = log2(mean(c(`Ratio H/L normalized 1`, `Ratio H/L normalized 2`), na.rm = T))) %>%
  mutate(`Ratio H/L normalized 1` = log2(`Ratio H/L normalized 1`), `Ratio H/L normalized 2` = log2(`Ratio H/L normalized 2`),
         Rep1_occ = log2(Rep1_occ), Rep2_occ = log2(Rep2_occ)) %>%
  ungroup() %>%
  mutate(SYMBOL = mapIds(org.Dm.eg.db,
                         keys=sub(";.*","", Proteins),
                         column="SYMBOL",
                         keytype="UNIPROT",
                         multiVals="first")) %>%
  filter(!is.na(Rep1_occ) | !is.na(Rep2_occ))

write_tsv(phospho_occupancy, file = "SuppTables/S3.txt")

correl <- round(cor(phospho_occupancy$Rep1_occ, phospho_occupancy$Rep2_occ, use = "complete.obs"), digits=4)

cairo_pdf("Figures/Figure-S3C.pdf")
ggplot(data = phospho_occupancy, aes(x = Rep1_occ, y = Rep2_occ))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_abline(slope = 1, intercept = 0, lty = "dotted")+
  theme_classic()+
  geom_point(data = phospho_occupancy[abs(phospho_occupancy$Rep1_occ) < 1 | abs(phospho_occupancy$Rep2_occ) < 1,],
             color = "black", alpha = 0.3)+
  geom_point(data = phospho_occupancy[abs(phospho_occupancy$Rep1_occ) >= 1 & abs(phospho_occupancy$Rep2_occ) >= 1,],
             color = "firebrick3", alpha = 1)+
  geom_text_repel(
    data = phospho_occupancy[abs(phospho_occupancy$Rep1_occ) >= 1 & 
                               abs(phospho_occupancy$Rep2_occ) >= 1,],
    aes(x = Rep1_occ, y = Rep2_occ, label = SYMBOL),
    size = 2)+
  coord_fixed(ratio = 1/1)+
  annotate(geom="text", x=-1, y=4, label= paste("R = ", correl), color="black")+
  scale_x_continuous(breaks = seq(-5,5,by = 1), name = "Replicate 1 (log2 Occupancy Holidic/Yeast food)") +
  scale_y_continuous(breaks = seq(-5,5,by = 1), name = "Replicate 2 (log2 Occupancy Holidic/Yeast food)")
dev.off()

## -- Figure S3D -- ##
### Prepare a ranked list
phospho_entrez <- UniProt.ws::select(up, sub(";.*", "", phospho_occupancy$Protein), "ENTREZ_GENE") %>%
  distinct(UNIPROTKB, .keep_all = T)

phospho_occupancy %>%
  rowwise() %>%
  mutate(mean =mean(c(Rep1_occ, Rep2_occ), na.rm = T)) %>%
  ungroup() %>%
  mutate(Proteins_simple = sub(";.*", "", phospho_occupancy$Protein)) %>%
  left_join(phospho_entrez, by = c("Proteins_simple" = "UNIPROTKB")) %>%
  dplyr::select(ENTREZ_GENE, mean)  %>%
  mutate(ENTREZ_GENE = as.numeric(as.character(ENTREZ_GENE))) -> phospho_gsea

phospho_gsea %>%
  mutate(mean = 2^mean) %>%
  arrange(mean) -> phospho_gsea_webgestalt
  
write.table(file = "Data/Figure-2/GSEA_phospho/Holidic_phospho_gsea_nonlog2.rnk", phospho_gsea_webgestalt, sep = "\t", col.names = F, row.names = F, quote = F)

cairo_pdf("Figures/Figure-S3D.pdf", 5, 5)
plot_gsea(entrez_gene_kegg_annot = "Libraries/dmelanogaster_pathway_KEGG_entrezgene.gmt",
          quant_data = "Data/Figure-2/GSEA_phospho/Holidic_phospho_gsea_KEGG/interestingID_mappingTable_wg_result1610578564.txt",
          gsea_result = "Data/Figure-2/GSEA_phospho/Holidic_phospho_gsea_KEGG/enrichment_results_wg_result1610578564.txt",
          pValue = 0.05)
dev.off()

## -- Figure S3E -- ##
which_element <- function(x) {
  protein_vector <- str_split(x["Proteins"], ";") %>% unlist()
  index <- which(protein_vector == x["Protein"])
  position_all <- str_split(x["Positions within proteins"], ";") %>% unlist()
  position <- as.numeric(as.character(position_all[index] %>% unlist()))
  return(position)
}

holidic_phospho <- read_tsv("Data/Figure-2/Phospho (STY)Sites.txt") %>%
  filter(`Localization prob`>=0.75) %>%
  filter(!grepl("REV__", Protein))
holidic_phospho[holidic_phospho == "NaN"] <- NA

holidic_phospho %>%
  mutate(Position_to_main = unlist(replace(apply(., 1, which_element), !sapply(apply(., 1, which_element), length),0))) %>%
  dplyr::select(Protein, Position_to_main, `Amino acid`, `Occupancy L 1`, `Occupancy H 1`, `Occupancy L 2`, `Occupancy H 2`) -> holidic_phospho_S3E

holidic_phospho_S3E[holidic_phospho_S3E == 0] <- NA
rowsums <- rowSums(is.na(holidic_phospho_S3E %>% dplyr::select(contains("Occupancy"))))

entrez_gene_kegg_annot = "Libraries/dmelanogaster_pathway_KEGG_entrezgene.gmt"
entrez_gene_kegg_annot <- read_lines(entrez_gene_kegg_annot)

entrez_gene_kegg_annot %>%
  map(str_split,pattern = "\\\t") %>% 
  map(unlist) %>% 
  map(1) -> kegg_ids_list

entrez_gene_kegg_annot %>%
  map2(.x = .,
       .y = kegg_ids_list,
       .f = kegg_ids_tab_to_long) %>% 
  map(as_tibble) %>% 
  map(separate_rows, value, sep = "\\\t") %>% 
  map2(.x = .,
       .y = kegg_ids_list,
       .f =  add_kegg_id) %>% 
  bind_rows() %>% 
  dplyr::rename(geneSet = kegg_id,
                Id = value) %>% 
  dplyr::select(geneSet, Id) -> entrez_gene_kegg_annot_long

gsea_result <- read_tsv("Data/Figure-2/GSEA_phospho/Holidic_phospho_gsea_KEGG/enrichment_results_wg_result1610578564.txt") %>%
  filter(pValue < 0.05)

cairo_pdf("Figures/Figure-S3E.pdf", 15, 5)

holidic_phospho_S3E %>%
  mutate(valid = rowsums) %>%
  filter(valid == 0) %>%
  gather(fraction_sample, occupancy, contains("Occupancy")) %>%
  mutate(fraction_sample = str_replace_all(fraction_sample, "Occupancy", "")) %>%
  mutate(fraction = str_replace_all(fraction_sample, " \\d", "")) %>%
  mutate(fraction = str_replace_all(fraction, " ", "")) %>%
  mutate(sample = as.numeric(str_replace_all(fraction_sample, "L|H ", ""))) %>%
  group_by(Protein, Position_to_main, fraction, `Amino acid`) %>%
  dplyr::summarise(mean_occ = mean(occupancy)) %>%
  ungroup() %>%
  mutate(Id = mapIds(org.Dm.eg.db,
                      keys=Protein,
                      column="ENTREZID",
                      keytype="UNIPROT",
                      multiVals="first")) %>%
  mutate(SYMBOL = mapIds(org.Dm.eg.db,
                     keys=Protein,
                     column="SYMBOL",
                     keytype="UNIPROT",
                     multiVals="first")) %>%
  left_join(entrez_gene_kegg_annot_long) %>%
  filter(geneSet %in% gsea_result$geneSet) %>%
  left_join(gsea_result) %>%
  mutate(fraction = fct_relevel(fraction, "L", "H")) %>%
  mutate(ID_plot = paste(SYMBOL, ", ", `Amino acid`, Position_to_main, sep = "")) -> holidic_phospho_S3E_tmp

holidic_phospho_S3E_tmp %>%
  group_by(ID_plot) %>%
  dplyr::summarise(n = n()/2) %>%
  left_join(holidic_phospho_S3E_tmp) %>%
  mutate(description = case_when(n > 1 ~ "Multiple", TRUE ~ description)) %>%
  mutate(description = fct_relevel(description, 
                                   "Protein processing in endoplasmic reticulum",
                                   "Ribosome",
                                   "Oxidative phosphorylation",
                                   "Glyoxylate and dicarboxylate metabolism",
                                   "Biosynthesis of amino acids",
                                   "Pyruvate metabolism",
                                   "Citrate cycle (TCA cycle)",
                                   "Multiple")) %>%
ggplot(aes(x = ID_plot, y = mean_occ, fill = fraction))+
  geom_bar(stat="identity", position=position_dodge())+
  facet_grid(.~ description, scales = "free", space = "free")+
  theme_bw()+
  scale_fill_manual(values = c("grey", "#CD2626"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2))

dev.off()

## -- Figure 2D -- ##
holidic_phospho_init <- read_tsv("Data/Figure-2/Phospho (STY)Sites.txt")
holidic_phospho_init[holidic_phospho_init == "NaN"] <- NA

print(paste("All detected phosphosites:", dim(read_tsv("Data/Figure-2/Phospho (STY)Sites.txt") %>%
                                                filter(!grepl("REV__", Protein)))[1]))
print(paste("Phosphosites p(loc) >= 0.75:", dim(read_tsv("Data/Figure-2/Phospho (STY)Sites.txt") %>%
                                                  filter(!grepl("REV__", Protein)) %>%
                                                  filter(`Localization prob`>=0.75))[1]))
print(paste("Occupancy quantified:", dim(read_tsv("Data/Figure-2/Phospho (STY)Sites.txt") %>%
                                           filter(!grepl("REV__", Protein)) %>%
                                           filter(`Localization prob`>=0.75) %>%
                                           mutate(Rep1_occ = log2(`Occupancy H 1`/`Occupancy L 1`), Rep2_occ = log2(`Occupancy H 2`/`Occupancy L 2`)) %>%
                                           filter(!is.na(Rep1_occ) | !is.na(Rep2_occ)))[1]))
print(paste("Unique phospho-protein:", dim(read_tsv("Data/Figure-2/Phospho (STY)Sites.txt") %>%
                                           filter(!grepl("REV__", Protein)) %>%
                                           filter(`Localization prob`>=0.75) %>%
                                           mutate(Rep1_occ = log2(`Occupancy H 1`/`Occupancy L 1`), Rep2_occ = log2(`Occupancy H 2`/`Occupancy L 2`)) %>%
                                           filter(!is.na(Rep1_occ) | !is.na(Rep2_occ)) %>%
                                           distinct(Protein))[1]))


pdf("Figures/Figure-2D.pdf", 5, 5)
phospho_occupancy %>%
  ggplot(aes(occ_mean))+
  geom_histogram(breaks = seq(-3.5,3.5, by = 0.1))+
  theme_classic()+
  scale_x_continuous(breaks = seq(-3,3,by = 1))
dev.off()

#########################
##-- Figure 3 and S4 --##
#########################

dm.mitocarta <- read_tsv("Data/Figure-3-4/Fly_Interactome_200515.txt") %>%
  filter(Mitochondrial.Process != "Translation (mt)")
oxphos_library <- read.table("Data/Figure-3-4/OXPHOS_components_Mm.Dm.txt", header = T, sep = "\t") # 

## Prepare the data frame
lrpprc_proteome <- read.table(file = "Data/Figure-3-4/proteinGroups.txt", header = T, sep = "\t") %>%
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+") %>%
  filter(Only.identified.by.site != "+") %>%
  dplyr::select(Protein.IDs, contains("Ratio.H.L.normalized.Exp_"), contains("Intensity.H.Exp"), contains("Intensity.L.Exp")) %>%
  dplyr::rename(Rep1 = "Ratio.H.L.normalized.Exp_01", Rep2 = "Ratio.H.L.normalized.Exp_02", 
         Rep3 = "Ratio.H.L.normalized.Exp_03", Rep4 = "Ratio.H.L.normalized.Exp_04") %>%
  mutate(Rep1 = log2(Rep1), Rep2 = -log2(Rep2), Rep3 = log2(Rep3), Rep4 = -log2(Rep4)) %>%
  column_to_rownames("Protein.IDs")

## Limma stats
fit <- lmFit(lrpprc_proteome %>% dplyr::select(contains("Rep")))
fit <- eBayes(fit)
limma.result <- topTable(fit,number = Inf,confint = TRUE, adjust.method = "fdr")

## Annotation
limma.result %>%
  rownames_to_column("UNIPROT_ID") %>%
  mutate(protein_simple = sub(";.*","", UNIPROT_ID)) %>%
  mutate(ENSEMBL = mapIds(org.Dm.eg.db, 
                         keys=as.character(protein_simple), 
                         column="ENSEMBL", 
                         keytype="UNIPROT",
                         multiVals="first")) %>%
  mutate(SYMBOL = mapIds(org.Dm.eg.db, 
                         keys=as.character(protein_simple), 
                         column="SYMBOL", 
                         keytype="UNIPROT",
                         multiVals="first")) %>%
  mutate(ENTREZ_ID = mapIds(org.Dm.eg.db, 
                         keys=as.character(protein_simple), 
                         column="ENTREZID", 
                         keytype="UNIPROT",
                         multiVals="first")) %>%
  left_join(dm.mitocarta, by = c("ENSEMBL" = "Gene.ID")) %>%
  mutate(mito = ENSEMBL %in% dm.mitocarta$Gene.ID) %>%
  left_join(lrpprc_proteome %>% rownames_to_column("UNIPROT_ID")) -> lrpprc_proteome

write.table(lrpprc_proteome, file = "SuppTables/Fig3.txt", sep = "\t", row.names = F, quote = F)

### Categorical t-test with fdr correction and boxplotting
cat.t.test <- function(data){
  cats <- unique(data$Mitochondrial.Process)
  cats <- cats[!is.na(cats)]
  p <- c()
  for(i in c(1:length(cats))){
    data.sub <- data[which(data$Mitochondrial.Process == cats[i]),]
    if(as.numeric(dim(data.sub)[1]) < 4){
      p.add <- NA
    } else {
      p.add <- p.adjust(unlist(t.test(data.sub$logFC)[3]), method = "fdr", n = length(data.sub$logFC))
    }
    p <- c(p,p.add)
  }
  stats <- data.frame(Category = cats, p.value = round(p,20))
  return(stats)
}

## -- Figure 3A -- ##
cat.t.test(lrpprc_proteome) -> t.test_lrpprc_proteome

write.table(t.test_lrpprc_proteome, file = "SuppTables/Fig3A_stats.txt", row.names = F, quote = F, sep = "\t")

cairo_pdf("Figures/Figure-3A.pdf", 10, 4)
lrpprc_proteome %>%
  drop_na(Mitochondrial.Process, logFC) %>%
ggplot(aes(x = fct_reorder(Mitochondrial.Process, logFC, .fun = median, .desc = T), y = logFC))+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_boxplot(outlier.shape = NA)+
  #coord_cartesian(ylim = c(-3,3))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_jitter(position=position_jitter(0.2), size = 0.4)+
  labs(title = "Mitochondrial functional categories in LRPPRC proteomics")+
  scale_y_continuous(breaks = seq(-5,3,by = 1))
dev.off()

## -- Figure 3B and S4D-- ##
lrpprc_proteome %>%
  drop_na(logFC) %>%
  left_join(oxphos_library %>% distinct(Dm, .keep_all = T), by = c("ENSEMBL" = "Dm")) %>%
  drop_na(Complex) %>%
  mutate(ID_complex = paste(Complex, ENSEMBL, SYMBOL.y, Mm, sep = "/")) %>%
  arrange(ID_complex) -> lrpprc_proteome_pheatmap_intermed

lrpprc_proteome_pheatmap_intermed %>%
  group_by(Complex) %>%
  dplyr::summarise(n = n()) -> gaps_row

gaps_row <- cumsum(gaps_row$n)

lrpprc_proteome_pheatmap_intermed %>%
  column_to_rownames("ID_complex") %>%
  dplyr::select(Rep1, Rep3, Rep2, Rep4) -> lrpprc_proteome_pheatmap

breakList <- seq(-2.1,2.1, by = 0.25)
colors <- colorpanel(length(breakList)-1, "deepskyblue", "black", "yellow")

cairo_pdf("Figures/Figure-3B.pdf", 5, 30)
pheatmap(lrpprc_proteome_pheatmap, color = colors, breaks = breakList, cluster_cols = F, cluster_rows = F,
         main = "LRPPRC fold changes per replicate", gaps_row = gaps_row)
dev.off()

## -- Figure S4A -- ##
holidic_overview <- holidic_proteome %>%
  mutate(ENSEMBL = mapIds(org.Dm.eg.db, 
                          keys=sub(";.*", "", Protein.IDs), 
                          column="ENSEMBL", 
                          keytype="UNIPROT",
                          multiVals="first")) %>%
  mutate(mitochondrial = ENSEMBL %in% dm.mitocarta$Gene.ID) %>%
  dplyr::select(contains("Intensity"), contains("mito")) %>%
  gather(Fraction, Intensity, 1:4) %>%
  mutate(type = "holidic")

lrpprc_overview <- lrpprc_proteome %>%
  dplyr::select(contains("Intensity"), "mito") %>%
  dplyr::rename("mitochondrial" = "mito") %>%
  gather(Fraction, Intensity, 1:8) %>%
  mutate(type = "lrpprc")

full_join(holidic_overview, lrpprc_overview) %>%
  filter(Intensity > 0) %>%
  drop_na(Intensity) %>%
  group_by(mitochondrial, type) %>%
  dplyr::summarise(median_int = median(Intensity)) %>%
  spread(mitochondrial, median_int) %>%
  mutate(ratio = `TRUE` / `FALSE`)

cairo_pdf("Figures/Figure-S4A_A.pdf")
full_join(holidic_overview, lrpprc_overview) %>%
  filter(Intensity > 0) %>%
  drop_na(Intensity) %>%
  ggplot(aes(x = type, y = log2(Intensity), fill = mitochondrial))+
  geom_boxplot()+
  coord_cartesian(ylim=c(15,45))+
  theme_classic()
dev.off()

bsf.levels <- lrpprc_proteome[lrpprc_proteome$SYMBOL %in% "bsf",]
bsf.levels <- bsf.levels[,grepl("Intensity.[LH].", colnames(bsf.levels))] %>%
  gather(Fraction, Intensity)%>%
  mutate(condition = c(rep(c("KD", "Control"),2), rep(c("Control", "KD"),2)))

bsf.levels %>%
  group_by(condition) %>%
  dplyr::summarise(median_bsf = median(Intensity)) %>%
  spread(condition, median_bsf) %>%
  mutate(ratio = KD / Control)

cairo_pdf("Figures/Figure-S4A-B.pdf")
ggplot(data = bsf.levels, aes(x = condition, y = log2(Intensity)))+
  geom_point()+
  theme_classic()+
  coord_cartesian(ylim=c(15,45))+
  theme_classic()
dev.off()

## -- Figure S4C -- ##
lrpprc_proteome %>%
  drop_na(Mitochondrial.Process) %>%
  dplyr::select("ENTREZ_ID", "logFC") -> lrpprc_proteome_gsea

write.table(lrpprc_proteome_gsea, file = "Data/Figure-3-4/GSEA/lrpprc_proteome_gsea.rnk", row.names = F, col.names = F, quote = F, sep = "\t")

cairo_pdf("Figures/Figure-S4C.pdf", 6, 5)
plot_gsea(entrez_gene_kegg_annot = "Libraries/dmelanogaster_pathway_KEGG_entrezgene.gmt",
          quant_data = "Data/Figure-3-4/GSEA/lrpprc_proteome_gsea_KEGG/interestingID_mappingTable_wg_result1611233580.txt",
          gsea_result = "Data/Figure-3-4/GSEA/lrpprc_proteome_gsea_KEGG/enrichment_results_wg_result1611233580.txt")
dev.off()

## -- Figure S4E -- ##
read.table("Data/Figure-3-4/LRPPRC_Mouse/LRPPRC_Mouse.txt", sep = "\t", header = T) %>%
  dplyr::rename(Mm = "Gene.names", Mm_FC020 = "logFC.2weeks", Mm_FC034 = "logFC.3.4weeks", Mm_FC050 = "logFC.5weeks", Mm_FC070 = "logFC.7weeks", Mm_FC100 = "logFC.10weeks") %>%
  dplyr::select(Mm, Mm_FC020, Mm_FC034, Mm_FC050, Mm_FC070, Mm_FC100) %>%
  right_join(oxphos_library) %>%
  left_join(lrpprc_proteome %>% dplyr:: select(ENSEMBL, logFC), by = c("Dm" = "ENSEMBL")) %>%
  dplyr::rename(Dm_FC = logFC) %>%
  drop_na(Dm_FC, Mm_FC100) %>%
  filter(Complex %in% c("1", "2", "3", "4", "5")) %>%
  gather(Type, logFC, contains("_FC")) -> kuehl_lrpprc

cairo_pdf("Figures/Figure-S4E.pdf", 7, 4)
ggplot(data = kuehl_lrpprc, aes(x = Complex, y = logFC, fill = Type))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_brewer(palette = "Purples")
dev.off()

################################
##-- Figure 4 and S5A,B,D --##
################################

read_tsv("Data/Figure-3-4/Phospho (STY)Sites.txt") %>%
  filter(`Localization prob`>=0.75) -> lrpprc_phospho

lrpprc_phospho[lrpprc_phospho =="NaN"] <- NA

print(paste("Mitochondrial phospho sites:",
            dim(lrpprc_phospho %>%
            mutate(ENSEMBL = mapIds(org.Dm.eg.db, 
                                  keys=Protein, 
                                  column="ENSEMBL", 
                                  keytype="UNIPROT",
                                  multiVals="first")) %>%
            left_join(dm.mitocarta, by = c("ENSEMBL" = "Gene.ID")) %>%
            drop_na(Mitochondrial.Process))[1]))

lrpprc_phospho %>%
  filter(is.na(`Potential contaminant`)) %>%
  filter(is.na(Reverse)) %>%
  filter(!grepl("REV", Protein)) %>%
  dplyr::select(Proteins, contains("Positions within"), contains("Occupancy L Exp_"), contains("Occupancy H Exp_")) %>%
  mutate(Rep1_occ = log2(`Occupancy H Exp_01`/`Occupancy L Exp_01`), Rep2_occ = log2(`Occupancy L Exp_02`/`Occupancy H Exp_02`),
         Rep3_occ = log2(`Occupancy H Exp_03`/`Occupancy L Exp_03`), Rep4_occ = log2(`Occupancy L Exp_04`/`Occupancy H Exp_04`)) %>%
  mutate(ID = paste(Proteins, `Positions within proteins`, sep = "::")) %>%
  column_to_rownames("ID") %>%
  dplyr::select(Rep1_occ, Rep2_occ, Rep3_occ, Rep4_occ) %>%
  filter(rowSums(is.na(.)) < 4) -> lrpprc_phospho

lrpprc_phospho[lrpprc_phospho =="NaN"] <- NA

fit <- lmFit(lrpprc_phospho)
fit <- eBayes(fit)
limma.result_lrpprc_phospho <- topTable(fit,number = Inf,confint = TRUE, adjust.method = "fdr") %>% rownames_to_column("ID") %>% mutate(Protein = sub("::.*", "", ID))

lrpprc_phospho %>%
  rownames_to_column("ID") %>%
  left_join(limma.result_lrpprc_phospho) %>%
  mutate(Proteins = str_replace_all(ID, "::.*", "")) %>%
  mutate(Pos_in_prot = str_replace_all(ID, ".*::", "")) %>%
  mutate(ENSEMBL = mapIds(org.Dm.eg.db, 
                          keys=sub(";.*", "", Proteins), 
                          column="ENSEMBL", 
                          keytype="UNIPROT",
                          multiVals="first")) %>%
  mutate(SYMBOL = mapIds(org.Dm.eg.db, 
                          keys=sub(";.*", "", Proteins), 
                          column="SYMBOL", 
                          keytype="UNIPROT",
                          multiVals="first")) %>%
  mutate(ENTREZ = mapIds(org.Dm.eg.db, 
                          keys=sub(";.*", "", Proteins), 
                          column="ENTREZID", 
                          keytype="UNIPROT",
                          multiVals="first")) %>%
  left_join(dm.mitocarta, by = c("ENSEMBL" = "Gene.ID")) -> lrpprc_phospho

write.table(lrpprc_phospho, sep = "\t", row.names = F, file = "SuppTables/Fig4.txt", quote = F)

## -- Figure S5B -- ##
cairo_pdf("Figures/Figure-S5B.pdf", 4, 4)
ggplot(data = lrpprc_phospho, aes(x = logFC, y = -log10(adj.P.Val), label = SYMBOL))+
  geom_point(data = filter(lrpprc_phospho, is.na(Mitochondrial.Process)), color = "grey70")+
  geom_point(data = filter(lrpprc_phospho, !is.na(Mitochondrial.Process), Mitochondrial.Process != "Oxidative Phosphorylation"), color = "black")+
  geom_text_repel(data = filter(lrpprc_phospho, adj.P.Val < 0.05, Mitochondrial.Process != "Oxidative Phosphorylation"), size = 2, max.overlaps = 20)+
  geom_text_repel(data = filter(lrpprc_phospho, adj.P.Val < 0.05, is.na(Mitochondrial.Process)), size = 2, color = "grey70", max.overlaps = 20)+
  theme_classic()+
  geom_hline(yintercept = -log10(0.05), lty = "dotted")
dev.off()

## -- Figure 4A -- ##
lrpprc_phospho_modspec <- read_tsv("Data/Figure-3-4/modificationSpecificPeptides.txt") %>%
  mutate(ENSEMBL = as.character(mapIds(org.Dm.eg.db, 
                          keys=sub(";.*", "", Proteins), 
                          column="ENSEMBL", 
                          keytype="UNIPROT",
                          multiVals="first"))) %>%
  mutate(SYMBOL = as.character(mapIds(org.Dm.eg.db, 
                         keys=sub(";.*", "", Proteins), 
                         column="SYMBOL", 
                         keytype="UNIPROT",
                         multiVals="first"))) %>%
  mutate(ENTREZ = as.character(mapIds(org.Dm.eg.db, 
                         keys=sub(";.*", "", Proteins), 
                         column="ENTREZID", 
                         keytype="UNIPROT",
                         multiVals="first"))) %>%
  left_join(dm.mitocarta, by = c("ENSEMBL" = "Gene.ID"))

lrpprc_phospho_ratios <- left_join(
  left_join(
    dm.mitocarta %>%
    group_by(Mitochondrial.Process) %>%
    dplyr::summarise(N.1protein = n()),
    
    lrpprc_phospho_modspec %>%
      filter(Modifications == "Unmodified") %>%
      group_by(Mitochondrial.Process) %>%
      dplyr::summarise(N.2Unmod = n())),
  
  lrpprc_phospho_modspec %>%
    filter(grepl("Phospho", Modifications)) %>%
    dplyr::group_by(Mitochondrial.Process) %>%
    dplyr::summarise(N.3Phospho = n())) 

scaling_factor <- lrpprc_phospho_ratios %>%
  gather(type, N, contains("N.")) %>%
  group_by(type) %>%
  dplyr::summarise(sum = sum(N, na.rm = T)) 

lrpprc_phospho_ratios %>%
  mutate(Ratio.1protein = N.1protein / scaling_factor[scaling_factor$type == "N.1protein",]$sum) %>%
  mutate(Ratio.2unmod = N.2Unmod / scaling_factor[scaling_factor$type == "N.2Unmod",]$sum) %>%
  mutate(Ratio.3phospho = N.3Phospho / scaling_factor[scaling_factor$type == "N.3Phospho",]$sum) %>%
  gather(Type_ratio, Ratio, contains("Ratio.")) -> lrpprc_phospho_ratios_plot

cairo_pdf(file = "Figures/Figure-4A.pdf", 10, 3)
ggplot(data=lrpprc_phospho_ratios_plot, aes(x=Mitochondrial.Process, y=Ratio, fill=Type_ratio)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## -- Figure 4B -- ##
phospho.detected_curated <- read.table(file = "Data/Figure-S4/OXPHOS.phospho.all_curated.txt", sep = "\t", header = T)

phospho.detected_curated$mutated.from <- sub("\\..*","",phospho.detected_curated$AA.Dm)
phospho.detected_curated$mutated.to <- sub("\\..*","",phospho.detected_curated$AA.Hs)
phospho.detected_curated$combo <- paste(phospho.detected_curated$mutated.from, ">", phospho.detected_curated$mutated.to, sep = "")

phospho.alluvial <- data.frame(table(phospho.detected_curated$combo))
phospho.alluvial$from <- sub(">.*", "", paste("f.", phospho.alluvial$Var1))
phospho.alluvial$to <- sub(".*>", "", phospho.alluvial$Var1)

is_alluvia_form(as.data.frame(phospho.alluvial), axes = 2:4, silent = TRUE)

cairo_pdf("Figures/Figure-4B.pdf")
ggplot(as.data.frame(phospho.alluvial),
       aes(y = Freq, axis1 = from, axis2 = to))+
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(aes(fill = from), stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum(width = 1/12, fill = "white", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  guides(fill = FALSE)
dev.off()

## -- Figure 4C -- ##
lrpprc_phospho %>%
  drop_na(Mitochondrial.Process) %>%
  left_join(phospho.detected_curated, by = c("ID" = "short"))  -> lrpprc_phospho_4C

lrpprc_phospho_4C[is.na(lrpprc_phospho_4C)] <- FALSE

cairo_pdf("Figures/Figure-4C.pdf", 5, 5)
ggplot()+
  geom_point(data = lrpprc_phospho_4C %>% filter(conserved != TRUE), aes(x = logFC, y = -log10(adj.P.Val)), color = "#BEBDBD")+
  geom_point(data = lrpprc_phospho_4C %>% filter(grepl("Oxidative Phosphorylation", Mitochondrial.Process) == TRUE), aes(x = logFC, y = -log10(adj.P.Val)), color = "#606262", size = 3)+
  geom_point(data = lrpprc_phospho_4C %>% filter(conserved == TRUE), aes(x = logFC, y = -log10(adj.P.Val)), color = "#D29878", size = 3)+
  geom_point(data = lrpprc_phospho_4C %>% filter(conserved.p == TRUE), aes(x = logFC, y = -log10(adj.P.Val)), color = "#CB6828", size = 3)+
  geom_point(data = lrpprc_phospho_4C %>% filter(conserved.p == TRUE), aes(x = logFC, y = -log10(adj.P.Val)), color = "#B13825", size = 3)+
  geom_text_repel(data = lrpprc_phospho_4C %>% filter(conserved == TRUE), aes(x = logFC, y = -log10(adj.P.Val), label = humanized))+
  geom_text_repel(data = lrpprc_phospho_4C %>% filter(conserved == FALSE, adj.P.Val < 0.05, grepl("Oxidative Phosphorylation", Mitochondrial.Process) == TRUE), aes(x = logFC, y = -log10(adj.P.Val), label = SYMBOL))+
  geom_hline(yintercept = -log10(0.05), lty = "dotted")+
  theme_classic()
dev.off()

## -- Figure S4B -- ##
lrpprc_coverage <- left_join(
  left_join(
    left_join(
      left_join(
        left_join(dm.mitocarta, dm.mitocarta %>% mutate(mito = Gene.ID %in% lrpprc_proteome$ENSEMBL)) %>%
          dplyr::rename(ENSEMBL = Gene.ID) %>%
          dplyr::select(ENSEMBL, mito) %>%
          group_by(mito) %>%
          dplyr::summarise(a.mitocarta = n()) %>%
          mutate(a.mitocarta = a.mitocarta/sum(a.mitocarta)),
        lrpprc_proteome %>%
          mutate(mito = ENSEMBL %in% dm.mitocarta$Gene.ID) %>%
          dplyr::select(ENSEMBL, mito) %>%
          group_by(mito) %>%
          dplyr::summarise(b.protein = n()) %>%
          mutate(b.protein = b.protein/sum(b.protein))),
      read_tsv("Data/Figure-3-4/Phospho (STY)Sites.txt") %>%
        filter(`Localization prob`>=0.75) %>%
        mutate(ENSEMBL = as.character(mapIds(org.Dm.eg.db, 
                                              keys=sub(";.*", "", Protein), 
                                              column="ENSEMBL", 
                                              keytype="UNIPROT",
                                              multiVals="first"))) %>%
        distinct(ENSEMBL) %>%
        mutate(mito = ENSEMBL %in% dm.mitocarta$Gene.ID)%>%
        group_by(mito) %>%
        dplyr::summarise(c.phosphoproteome = n()) %>%
        mutate(c.phosphoproteome = c.phosphoproteome/sum(c.phosphoproteome))),
    lrpprc_phospho_modspec %>%
      filter(Modifications == "Unmodified") %>%
      mutate(mito = ENSEMBL %in% dm.mitocarta$Gene.ID)%>%
      group_by(mito) %>%
      dplyr::summarise(d.unmodified_peptides = n()) %>%
      mutate(d.unmodified_peptides = d.unmodified_peptides/sum(d.unmodified_peptides))),
  lrpprc_phospho_modspec %>%
    filter(grepl("Phospho", Modifications)) %>%
    mutate(mito = ENSEMBL %in% dm.mitocarta$Gene.ID) %>%
    group_by(mito) %>%
    dplyr::summarise(e.modified_peptides = n()) %>%
    mutate(e.modified_peptides = e.modified_peptides/sum(e.modified_peptides))) %>%
  gather(type, ratio, 2:6)

cairo_pdf(file = "Figures/Figure-S4B.pdf", 5, 4)
ggplot(data=lrpprc_coverage, aes(x = type, y = ratio, fill = mito)) +
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## -- Figure S5A
phospho.detected_curated <- read.table(file = "Data/Figure-S4/OXPHOS.phospho.all_curated.txt", sep = "\t", header = T)
phospho.by.complex <- data.frame(complex = phospho.detected_curated$complex,
                                 conserved = phospho.detected_curated$conserved)

phospho.by.complex$combo <- paste(phospho.by.complex$complex, phospho.by.complex$conserved, sep =".")
phospho.by.complex <- data.frame(table(phospho.by.complex$combo))
phospho.by.complex$complex <- sub("\\..*","",phospho.by.complex$Var1)
phospho.by.complex$conserved <- sub(".*\\.","",phospho.by.complex$Var1)

cairo_pdf("Figures/Figure-S5A.pdf", 6, 4)
ggplot(data = phospho.by.complex, aes(x = complex, y = Freq, fill = conserved))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(breaks = seq(0,14, by = 2))
dev.off()

## -- Figure S5D -- ##
NDUFB10 <- lrpprc_phospho_modspec[lrpprc_phospho_modspec$Sequence=="YGDLGGYANAK",]
NDUFB10.int <- NDUFB10[,grepl("Ratio.H.L.normalized.", colnames(NDUFB10))]
NDUFB10.int[,c(2,4)] <- 1/NDUFB10.int[,c(2,4)]
NDUFB10.int <- data.frame(t(NDUFB10.int))
colnames(NDUFB10.int) <- c("unmod", "mod")

NDUFB10.int$Sample <- rownames(NDUFB10.int)
NDUFB10.int <- gather(NDUFB10.int, "Status", "HL", 1:2)

cairo_pdf("Figures/Figure-S5D.pdf",2,4)
ggplot(data = NDUFB10.int, aes(x = Status, y = log2(HL)))+
  geom_point()+
  theme_bw()+
  coord_cartesian(ylim=c(-2,0))
dev.off()
