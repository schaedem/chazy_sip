library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(RColorBrewer)
library(markdown)
library(ggstar)
library(ggplot2)
library(plotrix)
library(ggnewscale)
library(TDbook) 
library(phyloseq)
library(cluster)
library(dplyr)
library(tidyverse)
library(rcartocolor)
library(kableExtra)
library(ggtext)
library(ggsci)
library(cowplot)
library(ggpubr)


#normalized abundance values of incorporator ASVs with metadata and taxonomy; present in rarefied microcosm communities
incorp_norm <- readRDS("incorp_bulk_rare_abund.RDS") 
head(incorp_norm)

#identify dual-incorporating ASVs at each timepoint
dual_incorp <- incorp_norm %>%
  group_by(ASV, land, day, rep) %>%
  summarize(dual_incorp = sum(length(unique(labeled_substrate))),
            shared_incorp = sum(length(unique(land)))) %>%
  mutate(dual_incorp = ifelse(dual_incorp == "1", "0", 
                              ifelse(dual_incorp == "2", "1", dual_incorp))) %>%
  ungroup() %>%
  unique()

#merge back with main dataframe
incorp_norm_merge <- merge(incorp_norm, dual_incorp, by=c("ASV", "land", "day", "rep")) %>%
  unique() %>%
  group_by(ASV, land, day, rep) %>%
  mutate(labeled_substrate2 = ifelse(dual_incorp == 1, "Dual incorporator", labeled_substrate),
         c_sources = ifelse(dual_incorp ==1, 2, 1)) %>%
  ungroup()%>%
  mutate(labeled_substrate = ifelse(labeled_substrate == "13C", "Cellulose", "Xylose"), 
         land = ifelse(land == "NTH", "No till", "Plow till"))


##### Phyloseq object for incorporators ####
tree <- readRDS("tree_phy.RDS")
seqs <-  Biostrings::readDNAStringSet("dna-sequences.fasta")

meta <- read_csv("chazy.metadata.pico.csv") %>%
  dplyr::select(-c(Soil_Moisture, Soil_Temp, Sample_Date, Month, Year, Number, Fraction, 
                   PrimerF, PrimerR, Barcode, Rev.index, Fwd.index, Experiment, Density,
                   OM,  mean_pico_conc, adj_pico_conc)) %>% #get moisture, temp, pH data from 09-2014 metadata file
  filter(Full.sample %in% incorp_norm$sample) %>%
  column_to_rownames("Full.sample")

samp_data <- sample_data(meta)

otu_mat <- incorp_norm_merge %>%
  dplyr::select(sample, ASV, adj_abund_pap) %>%
  unique() %>%
  pivot_wider(names_from = sample, values_from = adj_abund_pap, values_fill = 0) %>%
  column_to_rownames("ASV") %>%
  as.matrix()

otu_mat <- otu_table(otu_mat, taxa_are_rows = T) 

taxa <- incorp_norm_merge %>%
  dplyr::select(ASV, Domain, Phylum, Class, Order, Family, Genus) %>%
  unique() %>% #337 unique incorporator ASVs (rarefied, present in bulk SIP microcosms)
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("ASV") %>%
  as.matrix()

incorp_taxa <- tax_table(taxa)
incorp_phylo <- phyloseq(otu_mat, incorp_taxa, samp_data, tree, seqs)


##### Identify abundant incorporators from each tillage regime ####
#group by variables to find most abundant taxa for plotting
ps.m = incorp_phylo %>%
  psmelt() %>%
  group_by(OTU, Land_Management, Day) %>%
  summarise(mAbund = mean(Abundance)) %>%
  ungroup() %>% 
  group_by(OTU, Land_Management) %>%
  summarise(smAbund = sum(mAbund)) 
head(ps.m)

#Top 35 most abundant taxa per soil type
topTax = ps.m %>% 
  group_by(Land_Management) %>%
  slice_max(order_by = smAbund, n = 35) %>%
  .$OTU %>%
  as.character()
length(unique(topTax)) #62 unique ASVs

#### Figure 6 #####

tax_df <- incorp_norm_merge2 %>%
  filter(ASV %in% topTax) %>%
  mutate(GenusLabel = ifelse(!is.na(Genus), paste(Genus),
                             ifelse(!is.na(Family), paste('Uncl. ', Family, sep = ""),
                                    ifelse(!is.na(Order), paste('Uncl. ', Order, sep = ""),
                                           ifelse(!is.na(Class), paste('Uncl. ', Class, sep = ""), paste("Uncl. ", Phylum, sep = "")))))) %>%
  mutate(GenusLabel = ifelse(GenusLabel == "Rhizobiaceae_putative  Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", 'Allo-Neo-Para-Rhizobium', GenusLabel),
         GenusLabel = ifelse(GenusLabel == "putative  Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", 'Allo-Neo-Para-Rhizobium', GenusLabel)) %>%
  mutate(GenusLabel = ifelse(GenusLabel == "unclassified  KD4-96", "Chloroflexi KD4-96", GenusLabel),
         GenusLabel = ifelse(GenusLabel == "putative  uncultured", paste("Uncl.", Family, sep=" "), GenusLabel),
         GenusLabel = ifelse(GenusLabel == "unclassified  uncultured", paste("Uncl.", Order, sep=" "), GenusLabel),
         GenusLabel = ifelse(GenusLabel == "uncultured", paste("Uncl.", Order, sep = " "), GenusLabel),
         GenusLabel = ifelse(GenusLabel == "Pir4 lineage", "Pirellulaceae Pir4", GenusLabel),
         GenusLabel = str_replace(GenusLabel, "^unclassified  (.*)", "Uncl. \\1"),
         GenusLabel = str_replace(GenusLabel, "putative  (.*)", "Uncl. \\1"),
         GenusLabel = ifelse(GenusLabel == "67-14", "Solirubrobacterales 67-14", GenusLabel),
         GenusLabel = ifelse(GenusLabel == "Gitt-GS-136", "Chloroflexi Gitt-GS-136", GenusLabel),
         GenusLabel = ifelse(GenusLabel == "MB-A2-108", "Actinobacteria MB-A2-108", GenusLabel)) %>%
  group_by(land, day, labeled_substrate2, GenusLabel, Phylum, ASV) %>%
  summarize(abund = mean(adj_abund_pap)) %>%
  ungroup() %>%
  mutate(Phylum = factor(Phylum, levels=c("Crenarchaeota", "Acidobacteriota", "Actinobacteriota", "Bacteroidota", "Chloroflexi", "Firmicutes", "Planctomycetota", "Proteobacteria")),
         ASV = fct_reorder(ASV, desc(Phylum))) 

tax_df_order <- tax_df[with(tax_df, order(desc(Phylum), desc(GenusLabel))), ]

#Color palette for highlighting by Phylum
phylum_pal <- RColorBrewer::brewer.pal(length(unique(tax_df_order$Phylum)), "Dark2")

text_cols <- tax_df_order %>%
  dplyr::select(ASV, Phylum) %>%
  unique() %>%
  mutate(phy_color = as.integer(Phylum),
         phy_color = phylum_pal[phy_color])
text_cols <- as.character(text_cols$phy_color)

#Top taxa plot
plot <- tax_df_order %>%
  ggplot(aes(x=as.factor(day), y=fct_inorder(ASV), size=abund, fill=labeled_substrate2)) +
  facet_wrap(~land) +
  scale_y_discrete(breaks=tax_df$ASV, labels=tax_df$GenusLabel) +
  theme_linedraw()  +
  xlab("Days since C addition") +
  ylab("") +
  scale_fill_manual(values = c("#3C5488FF",  "#F39B7FFF","#00A087FF"), name="Labeled substrate") +
  theme(#axis.text.y = element_text(face="italic", size=15),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=12),
    axis.title = element_text(size=15),
    strip.text = element_text(size=18),
    legend.text = element_text(size=12),
    legend.title = element_text(size=15)) +
  scale_size(range=c(3,15), name = "Mean normalized\nabundance") +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  annotate("rect", fill = phylum_pal[1], alpha = 0.1, xmin = -Inf, xmax = Inf, ymin = 61.5, ymax = Inf)  + #Crenarchaeota
  annotate("rect", fill = phylum_pal[2], alpha = 0.1, xmin = -Inf, xmax = Inf, ymin = 57.5, ymax = 61.5)  + #Acidobacteria
  annotate("rect", fill=phylum_pal[3], alpha=0.1, xmin=-Inf, xmax=Inf, ymin=46.5 , ymax=57.5 )+ #Actinobacteria
  annotate("rect", fill=phylum_pal[4], alpha=0.2, xmin=-Inf, xmax=Inf, ymin=45.5 , ymax=46.5 )+ #Bacteroidota
  annotate("rect", fill=phylum_pal[5], alpha=0.1, xmin=-Inf, xmax=Inf, ymin=41.5 , ymax=45.5 )+ #Chloroflexi
  annotate("rect", fill=phylum_pal[6], alpha=0.1, xmin=-Inf, xmax=Inf, ymin=38.5 , ymax=41.5 )+ #Firmicutes
  annotate("rect", fill=phylum_pal[7], alpha=0.2, xmin=-Inf, xmax=Inf, ymin=36.5 , ymax=38.5 )+ #Planctomycetota
  annotate("rect", fill=phylum_pal[8], alpha=0.1, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=36.5) + #Proteobacteria
  geom_point(shape=21, alpha=0.9)

#Color coded labels by taxa; combined with top taxa plot in illustrator software
labels <- ggplot(tax_df_order, aes(x=as.factor(day), y=fct_inorder(ASV), color=Phylum)) +
  scale_y_discrete(breaks=tax_df$ASV, labels=tax_df$GenusLabel) +
  theme_void()  +
  geom_text(aes(label = GenusLabel, color = Phylum), x = 0.5, size = 8, fontface="italic", hjust=1) +
  scale_color_manual(values=phylum_pal) +
  theme(legend.position = "none")

