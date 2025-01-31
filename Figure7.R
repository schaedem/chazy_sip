library(DESeq2)
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggtext)
library(doParallel) #to run multiple cores
library(ggsci)
library(ggtext)
library(ggpmisc)
library(emmeans)
library(VennDiagram)
library(cowplot)
library(gghalves)
library(ape)

incorp_norm <- readRDS("incorp_bulk_rare_abund.RDS") #normalized, rarefied abundance of incorporators in bulk microcosm soil

#### Format data to calculate max log2fold change (growth dynamics) ####
incorp_norm2 <- incorp_norm %>%
  dplyr::select(ASV, labeled_substrate, microcosm_substrate, land, day, rep, sample, adj_abund2, log_adj_abund, rep)

traits <- incorp_norm %>%
  dplyr::select(ASV, n16S, GC) %>%
  unique()

max_value <- incorp_norm2 %>%
  dplyr::rename(Replicate = rep) %>%
  group_by(ASV, land, Replicate, labeled_substrate, microcosm_substrate) %>%
  summarize(max_value = max(adj_abund2)) 

#baseline abundance in bulk microcosms in H2O only treatments
bulk.sip.phylo <- readRDS("/home/rstudio/files/pool5/pool5.RDS")
rare.bulk.sip.phylo <- rarefy_even_depth(bulk.sip.phylo, sample.size = 1567, rngseed = 1064) #losing one sample with 968 counts  

#read in DNA yield concentration datafor bulk microcosms; concentration unit is ng/ul
dat <- read_csv("/home/rstudio/files/chazy_sip_bulk_conc.csv") #112
map <- read_csv("/home/rstudio/files/chazy_sip_bulk_mapping.csv") #112

df <- merge(dat, map, by="Number") %>%
  dplyr::rename(sample = Sample)

h2o_rrn <- rare.bulk.sip.phylo@otu_table %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  pivot_longer(!ASV, names_to = "sample", values_to = "value") %>%
  filter(grepl("H2", sample)) %>%
  filter(grepl(".Rep", sample)) %>%
  filter(ASV %in% max_value$ASV) %>% #filter H2O baselines for incorporator ASVs with non-zero max abundance
  merge(df, by="sample") %>% #get DNA yield info
  merge(traits) %>%
  mutate(counts = 1567) %>%
  mutate(adj_conc = (conc*100)/0.5, #conc is in unit ng/ul, and we want to adjust it to ng/g soil; multiply by 100 bc that is the total extraction volume; then divide by 0.5 g (mass of soil extracted)
         rel_abund = value/counts,
         adj_abund = rel_abund * adj_conc / n16S)
head(h2o_rrn)

h2o_sum_rrn <- h2o_rrn %>%
  group_by(sample) %>%
  summarize(adj_rrn = n16S*counts*rel_abund,
            sum_rrn = sum(adj_rrn)) %>%
  ungroup() %>%
  dplyr::select(sample, sum_rrn) %>%
  unique()

#Obtain final normalized abundance by dividing by total estimated rrn in the sample
h2o_abund <- h2o_rrn %>%
  merge(h2o_sum_rrn) %>%
  mutate(adj_abund2 = adj_abund/sum_rrn,   
         log_adj_abund = ifelse(adj_abund2>0, log(adj_abund2), 0))
head(h2o_abund)

#finally, assign baseline normalized abundance
h2o_abund <- h2o_abund %>%
  dplyr::rename(land = Till) %>%
  group_by(ASV, land, Replicate) %>%
  summarize(baseline_value = mean(adj_abund2),
            baseline_count = mean(value)) %>% #baseline is calculated as the mean normalized abundance across reps
  ungroup() 

##### calculating log2fold change in NTH using DESeq2 ####
fold_change <- merge(max_value, h2o_abund, by=c("ASV", "land", "Replicate")) %>%
  merge(traits)

#NTH log-fold change
coldata_nth <- fold_change %>%
  filter(land == "NTH") %>%
  dplyr::select(ASV, land, Replicate, labeled_substrate, microcosm_substrate, max_value, baseline_value, baseline_count, n16S) %>%
  pivot_longer(cols = c("max_value", "baseline_value"), names_to = "condition", values_to = "value") %>%
  mutate(sample = paste(land, labeled_substrate,  sep="_"),
         sample = paste(sample, microcosm_substrate, sep = "_"),
         sample = paste(sample, Replicate, sep = "_"),
         sample = paste(sample, condition, sep="_")) %>%
  dplyr::select(condition, sample, labeled_substrate, microcosm_substrate, Replicate, land) %>%
  unique() %>%
  column_to_rownames("sample") %>%
  mutate(Replicate = as.factor(Replicate), condition = as.factor(condition))
#Baseline abundance
baseline_counts_nth <- fold_change %>%
  filter(land == "NTH") %>%
  dplyr::select(ASV, land, Replicate, labeled_substrate, microcosm_substrate, baseline_count) %>%
  mutate(condition = "baseline_value") %>%
  mutate(sample = paste(land, labeled_substrate,  sep="_"),
         sample = paste(sample, microcosm_substrate, sep = "_"),
         sample = paste(sample, Replicate, sep = "_"),
         sample = paste(sample, condition, sep="_")) %>%
  dplyr::select(ASV, sample, baseline_count) %>%
  pivot_wider(names_from = "sample", values_from="baseline_count", values_fill=0) 
#Max observed abundance
count_data_nth <- incorp_norm %>%
  filter(land == "NTH") %>%
  mutate(sample = paste(land, labeled_substrate,  sep="_"),
         sample = paste(sample, microcosm_substrate, sep = "_"),
         sample = paste(sample, rep, sep = "_"),
         sample = paste(sample, "max_value", sep="_")) %>%
  group_by(ASV, sample) %>%
  summarize(max_value = max(value)) %>%
  pivot_wider(names_from = sample, values_from = max_value, values_fill = 0) %>%
  merge(baseline_counts_nth, by="ASV") %>%
  column_to_rownames("ASV") 

#rearrange so that rows in coldata are in same order as columns in count_data
count_data_ordered_nth <- count_data_nth[ ,order(match(colnames(count_data_nth), rownames(coldata_nth)))]

count_data_scaled_nth <- count_data_ordered_nth + 1 %>%
  as.matrix()

#Run Deseq for NTH data
dds <- DESeqDataSetFromMatrix(countData = count_data_scaled_nth,
                              colData = coldata_nth,
                              design = ~ condition)
dds <-DESeq(dds)
res <- results(dds)
resultsNames(dds)
res_nth <- results(dds, name="condition_max_value_vs_baseline_value") %>% as.data.frame()

res_rrn_nth <- res_nth %>%
  rownames_to_column("ASV") %>%
  merge(traits) %>%
  mutate(land = "NTH") %>%
  mutate(log2FoldChange = ifelse(log2FoldChange < 0, 0, log2FoldChange)) %>%
  filter(ASV %in% nonzero_abund$ASV)

#### Plow till max log-fold change ####
coldata_pth <- fold_change %>%
  filter(land == "PTH") %>%
  dplyr::select(ASV, land, Replicate, labeled_substrate, microcosm_substrate, max_value, baseline_value, baseline_count, n16S) %>%
  pivot_longer(cols = c("max_value", "baseline_value"), names_to = "condition", values_to = "value") %>%
  mutate(sample = paste(land, labeled_substrate,  sep="_"),
         sample = paste(sample, microcosm_substrate, sep = "_"),
         sample = paste(sample, Replicate, sep = "_"),
         sample = paste(sample, condition, sep="_")) %>%
  dplyr::select(condition, sample, labeled_substrate, microcosm_substrate, Replicate, land) %>%
  unique() %>%
  column_to_rownames("sample") %>%
  mutate(Replicate = as.factor(Replicate), condition = as.factor(condition))
#baseline abundance PTH
baseline_counts_pth <- fold_change %>%
  filter(land == "PTH") %>%
  dplyr::select(ASV, land, Replicate, labeled_substrate, microcosm_substrate, baseline_count) %>%
  mutate(condition = "baseline_value") %>%
  mutate(sample = paste(land, labeled_substrate,  sep="_"),
         sample = paste(sample, microcosm_substrate, sep = "_"),
         sample = paste(sample, Replicate, sep = "_"),
         sample = paste(sample, condition, sep="_")) %>%
  dplyr::select(ASV, sample, baseline_count) %>%
  pivot_wider(names_from = "sample", values_from="baseline_count", values_fill=0) 
#max observed abundance PTH
count_data_pth <- incorp_norm %>%
  filter(land == "PTH") %>%
  mutate(sample = paste(land, labeled_substrate,  sep="_"),
         sample = paste(sample, microcosm_substrate, sep = "_"),
         sample = paste(sample, rep, sep = "_"),
         sample = paste(sample, "max_value", sep="_")) %>%
  group_by(ASV, sample) %>%
  summarize(max_value = max(value)) %>%
  pivot_wider(names_from = sample, values_from = max_value, values_fill = 0) %>%
  merge(baseline_counts_pth, by="ASV") %>%
  column_to_rownames("ASV") 

#rearrange so that rows in coldata are in same order as columns in count_data
count_data_ordered_pth <- count_data_pth[ ,order(match(colnames(count_data_pth), rownames(coldata_pth)))]

count_data_scaled_pth <- count_data_ordered_pth + 1 %>%
  as.matrix()
#Deseq to calculate max log-fold change in PTH
dds <- DESeqDataSetFromMatrix(countData = count_data_scaled_pth,
                              colData = coldata_pth,
                              design = ~ condition)
dds <-DESeq(dds)
res_pth <- results(dds, name="condition_max_value_vs_baseline_value") %>% as.data.frame()

res_rrn_pth <- res_pth %>%
  rownames_to_column("ASV") %>%
  merge(traits) %>%
  mutate(land = "PTH") %>%
  mutate(log2FoldChange = ifelse(log2FoldChange < 0, 0, log2FoldChange)) %>%
  filter(ASV %in% nonzero_abund$ASV)

####  Merge NTH and PTH log2foldchange results ####
meta <- incorp_norm %>%
  dplyr::select(ASV, labeled_substrate, land) %>%
  mutate(labeled_substrate = ifelse(labeled_substrate == "13X", "Xylose", "Cellulose"),
         land = ifelse(land == "NTH", "No till", "Plow till"))

res_rrn <- rbind(res_rrn_nth, res_rrn_pth) %>%
  mutate(land = ifelse(land == "NTH", "No till", "Plow till")) %>%
  dplyr::select(ASV, log2FoldChange, n16S, land) %>%
  merge(meta) %>%
  unique() %>%
  dplyr::select(ASV, log2FoldChange, n16S, land, labeled_substrate)%>%
  dplyr::rename(max_l2fc = log2FoldChange, 
                substrate = labeled_substrate)


#### Calculate incorporator latency ####

incorp <- readRDS("hr_sip_incorp_parsed.RDS") %>%
  filter(padj < 0.05) %>% #1210
  unique() %>%
  mutate(min_peak = ifelse(substrate == "13X", 2, 6),
         latency = as.numeric(day)-min_peak,
         latency = ifelse(latency<0,0,latency)) %>%
  dplyr::select(ASV, log2FoldChange, substrate, latency, day, land, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  dplyr::rename(amt_label = log2FoldChange) %>% #degree of labeling is equal to shift in bouyant density of DNA after labeling
  mutate(land = ifelse(land == "NTH", "No till", "Plow till"))

#test normality
shapiro.test(incorp$latency) #non-normal
shapiro.test(log(incorp$latency+1)) #non-normal


#### Combine life traits for each incorporator ASV ####
incorp_asv <- incorp %>%
  group_by(ASV, land, substrate) %>%
  summarize(latency = median(latency),
            amt_label = median(amt_label)) %>%
  mutate(substrate = ifelse(substrate == "13C", "Cellulose", "Xylose"))

full_traits <- merge(res_rrn, incorp_asv)

full_traits %>%
  group_by(substrate, land) %>%
  summarize(mean = mean(max_l2fc),
            median = median(max_l2fc))

nth_traits <- full_traits %>%
  filter(land == "No till")

pth_traits <- full_traits %>%
  filter(land == "Plow till")

#### Figure 7 ####

substrate_latency <- ggplot(full_traits, aes(x=land, y=latency, fill=land)) +
  geom_violin(bw=2) +
  facet_wrap(~substrate) +
  stat_summary(fun.data=mean_sd, geom="pointrange", color="black") +
  ggsci::scale_fill_npg() +
  theme_linedraw() +
  stat_compare_means(method="wilcox.test", label = "p.signif", size=10, vjust=2, label.x=1.5, hjust=0.5) +
  theme(strip.text = element_text(size=15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=15),
        legend.position = "none",
        plot.margin = unit(c(0, 0.25, 0, 0.25), 
                           "inches")) +
  ylab("Latency") +
  xlab("")

substrate_l2fc <- ggplot(full_traits, aes(x=land, y=max_l2fc, fill=land)) +
  geom_violin(bw=0.6) +
  facet_wrap(~substrate) +
  stat_summary(fun.data=mean_sd, geom="pointrange", color="black") +
  ggsci::scale_fill_npg() +
  theme_linedraw() +
  stat_compare_means(method="wilcox.test", label = "p.signif", size=10, vjust=2, label.x=1.5, hjust=0.5) +
  theme(strip.text = element_text(size=15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=15),
        legend.position = "none",
        plot.margin = unit(c(0, 0.25, 0, 0.25), 
                           "inches")) +
  ylab(expression('Maximum log'["2"]~"fold change")) +
  xlab("") 

substrate_assimilation <- ggplot(full_traits, aes(x=land, y=amt_label, fill=land))+
  geom_violin(bw=1) +
  facet_wrap(~substrate) +
  stat_summary(fun.data=mean_sd, geom="pointrange", color="black") +
  ggsci::scale_fill_npg() +
  theme_linedraw() +
  stat_compare_means(method="wilcox.test", label = "p.signif", size=10, vjust=2, label.x=1.5, hjust=0.5) +
  theme(strip.text = element_text(size=15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=15),
        legend.position = "none",
        plot.margin = unit(c(0, 0.25, 0, 0.25), 
                           "inches")) +
  ylab(expression('Degree of '^"13"*"C assimilation")) +
  xlab("") 

substrate_rrn <- ggplot(full_traits, aes(x=land, y=log(n16S), fill=land)) +
  geom_violin(bw=0.25) +
  facet_wrap(~substrate) +
  stat_summary(fun.data=mean_sd, geom="pointrange", color="black") +
  scale_fill_npg()+
  theme_linedraw() +
  stat_compare_means(method="wilcox.test", label = "p.signif", size=10, vjust=2, label.x=1.5, hjust=0.5) +
  theme(strip.text = element_text(size=15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=15),
        legend.position = "none",
        plot.margin = unit(c(0, 0.25, 0, 0.25), 
                           "inches")) +
  ylab(expression('Natural log predicted'~italic(rrn)))+
  xlab("")


Fig7 <- plot_grid(substrate_rrn, substrate_latency, substrate_assimilation, substrate_l2fc, labels="AUTO", ncol=1) 

#### wilcoxon tests for traits by substrate x tillage ####
cell_traits <- full_traits %>%
  filter(substrate == "Cellulose")

xyl_traits <- full_traits %>%
  filter(substrate == "Xylose")

head(cell_traits)

#max l2fc
cell_l2fc <- wilcox.test(max_l2fc~land, data=cell_traits) %>% broom::tidy() %>% mutate(substrate = "Cellulose", test = "max_l2fc")
xyl_l2fc <- wilcox.test(max_l2fc~land, data=xyl_traits) %>% broom::tidy() %>% mutate(substrate = "Xylose", test = "max_l2fc")

#rrn
cell_rrn <- wilcox.test(n16S~land, data=cell_traits) %>% broom::tidy() %>% mutate(substrate = "Cellulose", test = "rrn")
xyl_rrn <- wilcox.test(n16S~land, data=xyl_traits) %>% broom::tidy() %>% mutate(substrate = "Xylose", test = "rrn")

#13C assimilation
cell_lab <- wilcox.test(amt_label~land, data=cell_traits) %>% broom::tidy() %>% mutate(substrate = "Cellulose", test = "degree")
xyl_lab <- wilcox.test(amt_label~land, data=xyl_traits) %>% broom::tidy() %>% mutate(substrate = "Xylose", test = "degree")

#latency
cell_lat <- wilcox.test(latency~land, data=cell_traits) %>% broom::tidy() %>% mutate(substrate = "Cellulose", test = "latency")
xyl_lat <- wilcox.test(latency~land, data=xyl_traits) %>% broom::tidy() %>% mutate(substrate = "Xylose", test = "latency")

results <- rbind(cell_l2fc, xyl_l2fc, cell_rrn, xyl_rrn, cell_lab, xyl_lab, cell_lat, xyl_lat) %>%
  mutate(p_adjust = p.adjust(p.value),
         signif = ifelse(p_adjust < 0.05, "yes", "no"))

#median values to report

cell_traits %>%
  group_by(land) %>%
  summarize(rrn = median(n16S),
            latency = median(latency),
            amt_label = median(amt_label), 
            l2fc = median(max_l2fc))

xyl_traits %>%
  group_by(land) %>%
  summarize(rrn = median(n16S),
            latency = median(latency),
            amt_label = median(amt_label), 
            l2fc = median(max_l2fc))
