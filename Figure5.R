library(HTSSIP)
library(DESeq2)
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggtext)
library(doParallel) #to run multiple cores
library(ggsci)
library(ggtext)
library(ggpmisc)
library(VennDiagram)
library(eulerr)
library(ggpubr)
library(ggvenn)
library(ggVennDiagram)
library(cowplot)
library(ggplotify)

##### Figure 5 venn diagrams: detection vs. labeling ####
#we want to compare the number of labeled (p<0.05) to detected (all ASVs passing filtering) incorporators across tillage regimes
#read in raw DEseq results
l2fc <- read.csv("seq_l2fc.csv") %>%
  mutate(comparison = str_split_fixed(.id, "\\| ", 2)[,2],
         comparison = str_split_fixed(comparison, "\\==",4),
         substrate = comparison[,2],
         substrate = str_split_fixed(substrate, " & ", 2)[,1],
         substrate = str_sub(substrate, 2, -2),
         day = comparison[,3],
         day = str_split_fixed(day, " & ", 2)[,1],
         day = str_sub(day, 3, -2),
         land = comparison[,4],
         land = str_sub(land,2,-3)) %>%
  dplyr::rename(ASV = OTU) %>%
  dplyr::select(-comparison) 


#detection  
l2fc_df <- l2fc %>%
  unique() %>%
  group_by(ASV) %>%
  mutate(min_padj = min(padj)) %>%
  ungroup() %>%
  filter(min_padj < 0.05) %>%
  select(ASV, land) %>%
  unique() %>%
  mutate(value = as.logical("TRUE")) %>%
  tidyr::spread(key=land, value=value) %>%
  select(-ASV) %>%
  mutate(test = "Detection", 
         NTH = ifelse(is.na(NTH), as.logical(FALSE), NTH),
         PTH = ifelse(is.na(PTH), as.logical(FALSE), PTH)) %>%
  dplyr::rename(`No till` = NTH, `Plow till` = PTH)

pal <- c(PTH = "#4DBBD5FF", NTH = "#E64B35FF")
detection_plot = plot(euler(l2fc_df, by = list(test)),
                      legend = FALSE, labels = TRUE,
                      quantities = list(TRUE, fontsize=10), strip=TRUE, fills=pal) %>%
  ggplotify::as.ggplot()

detection_plot <- detection_plot + ggtitle("Detection") +
  theme(plot.title = element_text(hjust = 0.5, vjust=1, size=15, face="bold"))
detection_plot

#labeling
incorp_df = l2fc %>%
  filter(padj < 0.05) %>%
  select(ASV, land, substrate) %>%
  unique %>%
  mutate(value = as.logical("TRUE")) %>%
  tidyr::spread(key=land, value=value) %>%
  select(-ASV) %>%
  dplyr::rename("No till" = "NTH", "Plow till" = "PTH") %>%
  dplyr::mutate(substrate = ifelse(substrate == "13C", "Cellulose", "Xylose"), 
                `No till` = ifelse(is.na(`No till`), as.logical(FALSE), `No till`),
                `Plow till` = ifelse(is.na(`Plow till`), as.logical(FALSE), `Plow till`))
pal <- c("#E64B35FF", "#4DBBD5FF")

label_plot <- plot(euler(incorp_df, by = substrate),
                   legend = FALSE, labels = TRUE,
                   quantities = list(TRUE, fontsize=10), strip=TRUE, fills=pal)  %>%
  ggplotify::as.ggplot()

label_plot <- label_plot + ggtitle("Labeling") +
  theme(plot.title = element_text(hjust = 0.5,size=15, face="bold"),
        plot.margin = unit(c(1,1,1,1), "cm"))


venn_plots <- cowplot::plot_grid(detection_plot, label_plot, nrow=2)
