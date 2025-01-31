library(tidyverse)
library(vegan)
library(phyloseq)
library(stringr) 
library(ggsci)
library(ggtext)
library(lmerTest)
library(lme4)
library(RColorBrewer)
library(emmeans)
library(ggpubr)
library(cowplot)
library(gridtext)
library(Maaslin2)
library(ggtext)
setwd('/home/rstudio/files')

#phyloseq object containing rarefied taxa from bulk micocosm communities
#raw sequence data is availablee on NCBI SRA
rare.bulk.sip.phylo <-  readRDS("chazy_bulk_rare_phylo.RDS")

count_seqs <- function(pool){
  
  count_dat <- sample_sums(pool)
  names <- names(count_dat)
  values <- unname(count_dat)
  counts_1 <- cbind(names, values) %>% as.data.frame()
  colnames(counts_1) = c("sample", "counts") 
  
  counts_1 <- counts_1 %>%
    arrange(as.numeric(counts)) %>%
    mutate(counts = as.numeric(counts))
  counts_1$n <- seq.int(nrow(counts_1))
  
  plot <- ggplot(counts_1, aes(n, counts)) +
    geom_line()+
    theme_classic()
  
  return(list("counts"=counts_1, "plot"=plot))
  
}

#### Calculate relative abundance ####

rare.counts <- count_seqs(rare.bulk.sip.phylo@otu_table)

abund.rare <- rare.bulk.sip.phylo@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("Asv") %>%
  pivot_longer(-Asv, names_to = "sample") %>%
  merge(rare.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(label = str_split_fixed(sample, "\\.",4)[,1],
         land_management = str_split_fixed(sample, "\\.",4)[,2],
         day = str_split_fixed(sample, "\\.",4)[,3],
         rep = str_split_fixed(sample, "\\.",4)[,4],
         day = as.numeric(gsub('Day(\\d+$)', '\\1', day)),
         rep = as.numeric(gsub('Rep(\\d$)', '\\1', rep))) %>%
  mutate(day = ifelse(grepl("H2", sample), "H2O", day)) %>%
  #filter(!grepl("H2", sample)) %>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly
abund.rare %>% group_by(sample) %>% summarize(total = sum(rel_abund))

#get taxonomy
tax <- bulk.sip.phylo@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("Asv") 
head(tax)

#merge taxonomy with rarefied relative abundance
count_tax.rare <- merge(abund.rare, tax, by="Asv")
head(count_tax.rare)

#### Relative abundance plot - Fig 4A ####

# Get summary abundance of each taxa
taxa.meta.sum <- count_tax.rare %>%
  dplyr::select(-c(value, counts)) %>%
  pivot_longer(c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Asv"),
               names_to = "level", 
               values_to = "taxon") %>%
  mutate(land = recode_factor(land_management, "NTH" = "No till", "PTH" = "Plow till"),
         day = factor(day, levels=c("H2O", "1", "3", "7", "14", "30")))

#format phylum names
phylum_rel_abund <- taxa.meta.sum %>%
  filter(level == "Phylum") %>%
  group_by(taxon, land, day, rep, sample) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "*\\1*")) %>%
  ungroup()

#pool taxon <3% relative abundance into "other" category
taxon_pool <- phylum_rel_abund %>%
  group_by(land, day, taxon) %>%
  summarize(mean = mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) <3,
            mean=mean(mean), 
            .groups="drop") 

#arrange samples by relative abundance of proteobacteria
sample_order <- phylum_rel_abund %>%
  filter(taxon=="*Proteobacteria*") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order=1:nrow(.)) %>%
  dplyr::select(sample, order)

#arrange data frame to group samples by relative abundance within each taxon; 
#this will aid visual comparison between the barplots
join_taxon_rel_abund <- inner_join(phylum_rel_abund, taxon_pool, by="taxon") %>%
  mutate(taxon = ifelse(pool, "Other", taxon)) %>%
  group_by(sample, land, day, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = mean(mean), .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, rel_abund, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>%
  ungroup() %>%
  unique() %>%
  group_by(taxon) %>%
  arrange(desc(rel_abund), .by_group=TRUE) %>%
  inner_join(., sample_order, by="sample") %>% #arrange samples in decreasing order of Proteobacteria relative abundance
  mutate(sample = factor(sample),
         sample = fct_reorder(sample,order))

pal <- c("#A6CEE3","#1F78B4", "#B2DF8A", "#33A02C", "#D3D3D3","#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6" )

plow_till_rel_abund <- 
  join_taxon_rel_abund %>%
  filter(land == "Plow till") %>%
  ggplot(aes(x=sample, y=rel_abund, fill=taxon, group=day)) +
  geom_col(width=1) +
  facet_grid(~day, scale="free_x", space="free", switch="x") +
  #scale_fill_brewer(palette = "Paired") +
  scale_fill_manual(name=NULL, values=pal)+
  theme_classic() +
  theme(legend.text=element_markdown(size=12),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        title = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15, color="white"),
        plot.title = element_textbox_simple(
          size=20,
          padding = margin(5.5, 5.5, 5.5, 5.5),
          margin=margin(0,0,0,0),
          halign = 0.5,
          fill="black", color="white")) +
  scale_y_continuous(trans="reverse",expand=c(0,0), labels=c("100", "75", "50", "25", "0")) +
  ylab("Relative abundance (%)") +
  xlab("Time since substrate addition (days)") +
  labs(title = "Plow till", color="white")


notill_rel_abund <- 
  join_taxon_rel_abund %>%
  filter(land == "No till") %>%
  ggplot(aes(x=sample, y=rel_abund, fill=taxon, group=day)) +
  geom_col(width=1) +
  facet_grid(~day, scale="free_x", space="free", switch="x") +
  scale_fill_manual(name=NULL, values=pal)+
  theme_classic() +
  theme(legend.text=element_markdown(size=12),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        title = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15, color="white"),
        plot.title = element_textbox_simple(
          size=20,
          padding = margin(5.5, 5.5, 5.5, 5.5),
          margin=margin(0,0,0,0),
          halign = 0.5,
          fill="black", color="white"),
        legend.position = "bottom") +
  ylab("Relative abundance (%)") +
  xlab("Time since substrate addition (days)") +
  labs(title = "No till", color="white") +
  scale_y_continuous(trans="reverse",expand=c(0,0), labels=c("100", "75", "50", "25", "0"))

legend <- get_legend(notill_rel_abund)
final <- plot_grid(plow_till_rel_abund + theme(legend.position = "none"), 
                   notill_rel_abund + theme(legend.position="none"), 
                   nrow=1, labels = c("A",""))
final <- plot_grid(final, legend, nrow=2, rel_heights = c(2, 0.3))

#### Beta diversity plots  - Fig 4B-D ####
data <- abund.rare %>%
  dplyr::select(sample, Asv, value) %>%
  pivot_wider(names_from='Asv', values_from='value') %>%
  unique() %>%
  column_to_rownames("sample")

meta_df <- abund.rare %>%
  dplyr::select(sample, land_management, day, rep) %>%
  unique() %>%
  column_to_rownames("sample") %>%
  mutate(rep = factor(rep))

#### Bray-Curtis plots & permanova ####
bray_dist <- avgdist(data, sample=1567)
perm <- how(nperm=999)
set.seed(1234)
bray_perm <- adonis2(bray_dist ~ land_management*day, permutations=perm, data=meta_df)
bray_perm

#graph results
set.seed(12345)
bray_nmds <- metaMDS(bray_dist, distance=weighted_unifrac, k=3) %>%
  scores() %>%
  as_tibble(rownames="sample")%>%
  merge(meta) %>%
  mutate(day = factor(day, levels=c("H2O","1","3","7","14","30"))) %>%
  as_tibble() %>% unique() %>%
  mutate(day_land = paste( day,land_management, sep="_"),
         day_land = factor(day_land, levels=c("H2O_NTH", "1_NTH", "3_NTH", "7_NTH", "14_NTH", "30_NTH",
                                              "H2O_PTH", "1_PTH", "3_PTH", "7_PTH", "14_PTH", "30_PTH")),
         day_land = fct_relevel(day_land,"H2O_NTH", "1_NTH", "3_NTH", "7_NTH", "14_NTH", "30_NTH",
                                "H2O_PTH", "1_PTH", "3_PTH", "7_PTH", "14_PTH", "30_PTH")) %>%
  mutate(land_management = ifelse(land_management == "NTH", "No till", "Plow till"))

vec_df <- bray_nmds %>%
  group_by(land_management, day) %>%
  summarize(centroid_x = mean(NMDS1),
            centroid_y = mean(NMDS2)) 

pal2 <- c ("#F7C0B9","#C2E7F0","#EA6856","#6AC6DB","#F3A398","#A5DCE9",
           "#E64B35","#4DBBD5","#EE8577","#87D1E2","#FCDEDA","#E0F3F7","#E64B35", "#4DBBD5")
greyscale <- c("#F3F3F3", "#DEDEDE", "#CCCCCC", "#999999", "#666666", "#222222")
#dummy plot for grey scale
plot <- ggplot(bray_nmds, aes(x=NMDS1, y=NMDS1, fill=day, color=day)) +
  geom_point(shape=15, size=10) +
  scale_color_manual(values=greyscale) +
  labs(color="Day", shape="Day", label="Day", size="Day", fill="Day") +
  theme_linedraw() +
  guides(color=guide_legend(title="Day", nrow=1)) +
  theme(legend.position = "bottom",  legend.title = element_text(size=15),
        legend.text = element_text(size=12),
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,0,0),
  legend.justification = "center")
color_leg <- get_legend(plot)

#Bray-Curtis biplot
bray_nmds_plot <- ggplot(bray_nmds, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=day_land, shape=land_management), size=3)+
  scale_color_manual(values=pal2, name = "Tillage") +
  guides(color="none") +
  theme_linedraw() +
  theme(
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    axis.text = element_text(size=12),
    legend.title = element_text(size=15),
    legend.text = element_text(size=12),
    legend.position = "bottom",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(0,0,0,0),
    legend.justification = "center"
  ) +
  stat_ellipse(aes(linetype=land_management, color=land_management)) +
  guides(linetype = guide_legend(override.aes=list(color=c("#E64B35", "#4DBBD5"))), nrow=1) +
  labs(linetype = "Tillage", color="Tillage", label="Tillage", shape= "Tillage") +
  #NTH centroids
  geom_segment(aes(x=0.519, y=-0.183, xend = -0.0230, yend = 0.112), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  geom_segment(aes(x=-0.0230, y=0.112, xend=0.016, yend=0.0483), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  geom_segment(aes(x=0.016, y=0.0483, xend=0.0265, yend=0.0254), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  geom_segment(aes(x=0.0265, y=0.0254, xend=0.0626, yend=0.0184), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  geom_segment(aes(x=0.0626, y=0.0184, xend=0.132, yend=0.0355), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  #PTH centroids
  geom_segment(aes(x=0.483, y=-0.189, xend =-0.0938, yend=0.0196), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_segment(aes(x=-0.0938, y=0.0196, xend = -0.126, yend=0.00169), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_segment(aes(x=-0.126, y=0.00169, xend = -0.101, yend=-0.0439), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_segment(aes(x=-0.101, y=-0.0439, xend = -0.131, yend=-0.0222), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_segment(aes(x=-0.131, y=-0.0222, xend = -0.121, yend=-0.0186), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_richtext(aes(x=0.38, y=-0.33), label=c("**Tillage** *p=0.001*"), fill=NA, label.color=NA, col="black", hjust=0) +
  geom_richtext(aes(x=0.38, y=-0.38, label="**Day** *p=0.001*"), fill=NA, label.color=NA, hjust=0) +
  geom_richtext(aes(x=0.38, y=-0.42, label="Tillage * Day *p=0.172*"), fill=NA, label.color=NA, hjust=0)
bray_nmds_plot 

till_leg <- get_legend(bray_nmds_plot)

#### Weighted unifrac ####
weighted_unifrac <- UniFrac(rare.bulk.sip.phylo, weighted=TRUE)
set.seed(1234)
w_unifrac <- weighted_unifrac %>%
  as.matrix() 
w_unifrac_perm <- adonis2(w_unifrac ~ land_management*day, permutations=perm, data=meta_df)
w_unifrac_perm

#graph results
set.seed(12345)
w_nmds <- metaMDS(weighted_unifrac, distance=weighted_unifrac, k=3) %>%
  scores() %>%
  as_tibble(rownames="sample")%>%
  merge(meta) %>%
  mutate(day = factor(day, levels=c("H2O", "1", "3", "7", "14", "30"))) %>%
  as_tibble() %>% unique() %>%
  mutate(day_land = paste( day,land_management, sep="_"),
         day_land = factor(day_land, levels=c("H2O_NTH", "1_NTH", "3_NTH", "7_NTH", "14_NTH", "30_NTH",
                                              "H2O_PTH", "1_PTH", "3_PTH", "7_PTH", "14_PTH", "30_PTH")),
         day_land = fct_relevel(day_land,"H2O_NTH", "1_NTH", "3_NTH", "7_NTH", "14_NTH", "30_NTH",
                                "H2O_PTH", "1_PTH", "3_PTH", "7_PTH", "14_PTH", "30_PTH")) %>%
  mutate(land_management = ifelse(land_management == "NTH", "No till", "Plow till"))

vec_df <- w_nmds %>%
  group_by(land_management, day) %>%
  summarize(centroid_x = mean(NMDS1),
            centroid_y = mean(NMDS2)) 

#weighted unifrac biplot
w_nmds_plot <- ggplot(w_nmds, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=day_land, shape=land_management), size=3)+
  scale_color_manual(values=pal2) +
  theme_linedraw() +
  stat_ellipse(aes(linetype=land_management, color=land_management)) +
  theme(
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    axis.text = element_text(size=12),
    legend.position = "none"
  ) +
  #NTH centroids
  geom_segment(aes(x=0.367, y=0.0365, xend = -0.0824, yend = 0.0437), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  geom_segment(aes(x=-0.0824, y=0.0437, xend=-0.0133, yend=0.0077), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  geom_segment(aes(x=-0.0133, y=0.0077, xend=0.0267, yend=0.00363), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  geom_segment(aes(x=0.0267, y=0.00363, xend=0.0671, yend=-0.0289), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  geom_segment(aes(x=0.0671, y=-0.0289, xend=0.105, yend=-0.0379), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  #PTH centroids
  geom_segment(aes(x=0.385, y=-0.0188, xend =-0.128, yend=0.0307), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_segment(aes(x=-0.128, y=0.0307, xend = -0.119, yend=0.00199), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_segment(aes(x=-0.119, y=0.00199, xend = -0.0586, yend=0.00279), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_segment(aes(x=-0.0586, y=0.00279, xend = -0.0426, yend=-0.0140), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_segment(aes(x=-0.0426, y=-0.014, xend = -0.0471, yend=-0.00993), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_richtext(aes(x=0.26, y=-0.3), label=c("**Tillage** *p=0.001*"), fill=NA, label.color=NA, col="black", hjust=0) +
  geom_richtext(aes(x=0.26, y=-0.34, label="**Day** *p=0.001*"), fill=NA, label.color=NA, hjust=0) +
  geom_richtext(aes(x=0.26, y=-0.38, label="Tillage * Day *p=0.188*"), fill=NA, label.color=NA, hjust=0)
w_nmds_plot 

##### Unweighted unifrac ####
unweighted_unifrac <- UniFrac(rare.bulk.sip.phylo, weighted=FALSE)
set.seed(1234)
uw_unifrac <- unweighted_unifrac %>%
  as.matrix()
uw_unifrac_perm <- adonis2(uw_unifrac ~ land_management*day, permutations=perm, data=meta_df)
uw_unifrac_perm

#graph results
set.seed(12345)
uw_nmds <- metaMDS(unweighted_unifrac, distance = unweighted_unifrac, k=4) %>% 
  scores() %>%
  as_tibble(rownames = "sample") %>%
  merge(meta) %>%
  mutate(day = factor(day, levels=c("H2O", "1", "3", "7", "14", "30"))) %>%
  as_tibble() %>% unique() %>%
  mutate(day_land = paste( day,land_management, sep="_"),
         day_land = factor(day_land, levels=c("H2O_NTH", "1_NTH", "3_NTH", "7_NTH", "14_NTH", "30_NTH",
                                              "H2O_PTH", "1_PTH", "3_PTH", "7_PTH", "14_PTH", "30_PTH")),
         day_land = fct_relevel(day_land,"H2O_NTH", "1_NTH", "3_NTH", "7_NTH", "14_NTH", "30_NTH",
                                "H2O_PTH", "1_PTH", "3_PTH", "7_PTH", "14_PTH", "30_PTH")) %>%
  mutate(land_management = ifelse(land_management == "NTH", "No till", "Plow till"))

vec_df <- uw_nmds %>%
  group_by(land_management, day) %>%
  summarize(centroid_x = mean(NMDS1),
            centroid_y = mean(NMDS2)) 

#unweighted unifrac biplot
uw_nmds_plot <- ggplot(uw_nmds, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=day_land, shape=land_management), size=3)+
  scale_color_manual(values = pal2) +
  theme_linedraw() +
  stat_ellipse(aes(linetype=land_management, color=land_management)) +
  theme(
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    plot.title = element_text(size=20),
    axis.text = element_text(size=12),
    legend.position = "none"
  ) +
  #NTH centroids
  geom_segment(aes(x=-0.256, y=-0.152, xend = 0.0726, yend = 0.182), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  geom_segment(aes(x=0.0726, y=0.182, xend=-0.053, yend=0.0714), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  geom_segment(aes(x=-0.0530, y=0.0714, xend=-0.0614, yend=0.06), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  geom_segment(aes(x=-0.0614, y=0.06, xend=-0.117, yend=0.0369), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  geom_segment(aes(x=-0.117, y=0.0369, xend=-0.236, yend=-0.0679), arrow = arrow(length=unit(0.2, 'cm')), color="#E64B35FF", size=1.2) +
  #PTH centroids
  geom_segment(aes(x=-0.284, y=-0.204, xend = 0.196, yend=0.126), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_segment(aes(x=0.196, y=0.126, xend = 0.231, yend=0.0469), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_segment(aes(x=0.231, y=0.0469, xend = 0.0957, yend=-0.0477), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_segment(aes(x=0.0957, y=-0.0477, xend = 0.0634, yend=-0.0915), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_segment(aes(x=0.0634, y=-0.0915, xend = -0.0291, yend=-0.184), arrow = arrow(length=unit(0.2, 'cm')), color="#4DBBD5FF", size=1.2) +
  geom_richtext(aes(x=0.22, y=-0.3), label=c("**Tillage** *p=0.001*"), fill=NA, label.color=NA, col="black", hjust=0) +
  geom_richtext(aes(x=0.2, y=-0.35, label="**Day** *p=0.001*"), fill=NA, label.color=NA, hjust=0) +
  geom_richtext(aes(x=0.2, y=-0.4, label="**Tillage * Day** *p=0.045*"), fill=NA, label.color=NA, hjust=0)
uw_nmds_plot  

#arrange beta-diversity plots
legend2 <- plot_grid( till_leg,color_leg, nrow=1)
legend2
ggsave("nmds_leg.png", width=10, height=2, units="in")

plots <- plot_grid(bray_nmds_plot + theme(legend.position = "none"),
                   uw_nmds_plot, w_nmds_plot, labels = c("B", "C", "D"), nrow=1)
#add legend
plots_leg <- plot_grid(plots, legend2, nrow=2, rel_heights= c(4,0.3))

###combine alpha and beta diversity plots for Fig 4 ####

plot_grid(final, plots_leg, nrow=2, align = "hv", axis="bt", rel_heights=c(2,1))

##### Beta dispersion ####
#Are no-till communities more dispersed than plow-till?
#Bray-Curtis dispersion
meta <- meta %>% unique() %>% column_to_rownames("sample")
bd_bray <- betadisper(bray_dist, meta$land_management)
permutest(bd_bray) #p=0.001
anova(bd_bray) #p<0.00001

#Weighed Unifrac distance
w_unifrac <- as.dist(w_unifrac)
bd_w <- betadisper(w_unifrac, meta$land_management)
permutest(bd_w) #p=0.008
anova(bd_w) #p=0.01

#Unweighted Unifrac distance
uw_unifrac <- as.dist(uw_unifrac)
bd_uw <- betadisper(uw_unifrac, meta$land_management)
permutest(bd_uw) #p=0.09
anova(bd_uw) #p=0.1
