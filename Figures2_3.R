library(tidyverse)
library(lmerTest)
library(lme4)
library(cowplot)
library(ggsci) #for color palettes

setwd("/home/rstudio/files/gas_data/final_data")
data <- read_csv("c_min_df.csv")

#### Analysis 1: Do carbon mineralization rates and cumulative C mineralized differ between tillage? (12C + 13C) ####
#treat 12C and 13C as replicates since they were treated the same except for the isotope label
df <- data %>%
  filter(co2 == "TIC") %>%#interested in total carbon (12C + 13C) for this analysis
  filter(!grepl("H2O", sample)) %>%
  mutate(Treatment = factor(Treatment),
         Rep = factor(Rep),
         sample_day = factor(sample_day)) #land management and sample day are significant

#Posthoc pairwise comparisons 

#mixed model comparing total C mineralized over time
mod0 <- lmer(mg_cum_C ~ Land_Management*sample_day + (1|Rep), data=df)
anova(mod0)

#Posthoc means separation
posthoc.mod0 = lsmeans(mod0, pairwise~Land_Management*sample_day, adjust="bonferroni")

cum_means.df = data.frame(posthoc.mod0$lsmeans) %>%
  mutate(sample_day = as.numeric(as.character(sample_day)))

cum_pairwise.df = data.frame(posthoc.mod0$contrasts) %>%
  tidyr::separate(contrast, into=c("tillage_1", "day_1", "tillage_2", "day_2"), remove=FALSE) %>%
  filter(tillage_1 != tillage_2, day_1 == day_2) %>%
  mutate(day = as.numeric(as.character(day_1))) %>%
  dplyr::select(-day_1, -day_2) %>%
  mutate(contrast = paste(tillage_1, tillage_2, sep="-")) %>% #no significant differences at any timepoints
  dplyr::rename(sample_day = day)

#run a mixed model comparing tillage mineralization rates over time
mod1 <- lmer(C_rate_mg ~ Land_Management*sample_day  + (1|Rep), data=df)
summary(mod1)
anova(mod1)

# Run a posthoc pairwise comparisons and get the means
posthoc.mod1 = lsmeans(mod1, pairwise~Land_Management*sample_day, adjust="bonferroni")
summary(df$C_rate_mg)

rate_means.df = data.frame(posthoc.mod1$lsmeans) %>%
  mutate(sample_day = as.numeric(as.character(sample_day)))

rate_pairwise.df = data.frame(posthoc.mod1$contrasts) %>%
  tidyr::separate(contrast, into=c("tillage_1", "day_1", "tillage_2", "day_2"), remove=FALSE) %>%
  filter(tillage_1 != tillage_2, day_1 == day_2) %>%
  mutate(day = as.numeric(as.character(day_1))) %>%
  dplyr::select(-day_1, -day_2) %>%
  mutate(contrast = paste(tillage_1, tillage_2, sep="-")) %>% #significant contrasts between tillage only at days 1,2
  dplyr::rename(sample_day = day)

plot_rate <- filter(rate_pairwise.df, p.value<0.1) %>%
  mutate(label = ifelse(p.value < 0.05, "*", ".")) %>%
  mutate(mean_rate = 0.065, Land_Management = NA) %>%
  dplyr::select(mean_rate, Land_Management, sample_day, label)

#plot results
mod1_df <- df %>%
  group_by(Land_Management, sample_day) %>%
  summarize(mean_rate = mean(C_rate_mg),
            sd_rate = sd(C_rate_mg), 
            mean_cum = mean(mg_cum_C),
            sd_cum = sd(mg_cum_C)) %>%
  ungroup()

mod1_plot <- ggplot(mod1_df, aes(x=sample_day, y=mean_rate, color=Land_Management, group=Land_Management)) +
  geom_point(position=position_dodge(0.2)) +
  geom_line(position=position_dodge(0.2)) +
  theme_linedraw() +
  geom_errorbar(aes(ymin=mean_rate-sd_rate, ymax=mean_rate+sd_rate), width=0.2, position=position_dodge(0.2)) +
  scale_color_npg(labels=c("No till", "Plow till"), name="Tillage") +
  xlab("Time since substrate addition (days)") +
  labs(y=expression("C mineralization rate (mg/hr)")) +
  theme(legend.position = "top",
        strip.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12)) +
  geom_text(data=plot_rate, aes(x=sample_day+1, y=mean_rate, label=label), size=6, inherit.aes = FALSE) #for some reason, wanted to shift labels over; correcting for that by adding 1

mod1_plot



#total cumulative C plot
cum_C_plot <- ggplot(mod1_df, aes(x=sample_day,y=mean_cum, color=Land_Management, group=Land_Management)) +
  geom_point(position = position_dodge(0.5)) +
  theme_linedraw() +
  geom_line(position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin=mean_cum-sd_cum, ymax=mean_cum+sd_cum),position = position_dodge(0.5), width=0.2) +
  scale_color_npg(labels=c("No till", "Plow till"), name="Tillage") +
  xlab("Time since substrate addition (days)") +
  labs(y=expression("Cumulative C mineralized (mg)")) +
  theme(strip.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))

### H2O model for overall C mineralization rates
df2 <- data %>%
  filter(co2 == "TIC") %>%
  filter(grepl("H2O", sample)) %>%
  mutate(Treatment = factor(Treatment),
         Rep = factor(Rep),
         sample_day = factor(sample_day)) 

#run a mixed model comparing tillage mineralization rates over time
mod2 <- lmer(C_rate_mg ~ Land_Management*sample_day  + (1|Rep), data=df2)
anova(mod2) #only sample day is significant; no differenes by tillage

mod3 <- lmer(mg_cum_C~ Land_Management*sample_day  + (1|Rep), data=df2)
anova(mod2) #only sample day is significant; no differences by tillage


#plot H2O controls
mod2_df <- df2 %>%
  group_by(Land_Management, sample_day) %>%
  summarize(mean_rate = mean(C_rate_mg),
            sd_rate = sd(C_rate_mg),
            mean_cum = mean(mg_cum_C),
            sd_cum = sd(mg_cum_C)) %>%
  mutate(Land_Management= paste("H2O control", Land_Management, sep="-"))

both_mods <- rbind(mod1_df, mod2_df) %>%
  mutate(Land_Management = factor(Land_Management))

palette <- c("#a6a6a6", "#a6a6a6", ggsci::pal_npg()(2))

levels(both_mods$Land_Management)

#H2O control cumulative C mineralized plot
cum_C_plot <- ggplot(both_mods, aes(x=sample_day, y=mean_cum, color = Land_Management, group = Land_Management)) +
  geom_point() +
  geom_line() +
  theme_linedraw() +
  geom_errorbar(aes(ymin=mean_cum-sd_cum, ymax=mean_cum+sd_cum), width=0.2) +
  scale_color_manual(values = palette, labels=c( "H2O Control (no till)", "H2O Control (plow till)", "No till", "Plow till"), name="") +
  xlab("Time since substrate addition (days)") +
  labs(y=expression("Cumulative C mineralized (mg)")) +
  theme(legend.position = "top",
        strip.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))
cum_C_plot

#H2O control cumulative mineralization rate plot - Supplemental Figure
mods1_2_plot <-ggplot(both_mods, aes(x=sample_day, y=mean_rate, color=Land_Management, group=Land_Management)) +
  geom_point() +
  geom_line() +
  theme_linedraw() +
  geom_errorbar(aes(ymin=mean_rate-sd_rate, ymax=mean_rate+sd_rate), width=0.2) +
  scale_color_manual(values = palette, labels=c( "H2O Control (no till)", "H2O Control (plow till)", "No till", "Plow till"), name="") +
  xlab("Time since substrate addition (days)") +
  labs(y=expression("C mineralization rate (mg/hr)")) +
  theme(legend.position = "top",
        strip.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12)) +
  geom_text(data=plot_rate, aes(x=sample_day+1, y=mean_rate, label=label), size=6, inherit.aes = FALSE) #for some reason, wanted to shift labels over; correcting for that by adding 1
mods1_2_plot

plot_grid(cum_C_plot + theme(axis.title.x = element_blank()), mods1_2_plot + theme(legend.position = "none"), 
          align= "hv", axis = "bt", ncol=1, labels = "AUTO")


#### Analysis 2: Do 13C mineralization rates differ between tillage for each substrate? ####
#Perform analysis individually for each substrate
#Cellulose
df3 <- data %>%
  filter(!grepl("H2O", sample)) %>%
  filter(co2 == "45") %>%
  filter(Treatment == "13C") %>%
  mutate(Treatment = factor(Treatment),
         Rep = factor(Rep),
         sample_day = factor(sample_day)) 

mod3 <- lmer(C_rate_mg ~ Land_Management*sample_day  + (1|Rep), data=df3)
anova(mod3) #land management, day, and their interaction is significant

# Run a posthoc pairwise comparisons and get the means
#lsmeans(treat.model, pairwise~ecosystem*day, adjust="tukey")
posthoc.mod3 = lsmeans(mod3, pairwise~Land_Management*sample_day, adjust=NULL)

rate_C13.df = data.frame(posthoc.mod3$lsmeans) %>%
  mutate(sample_day = as.numeric(as.character(sample_day)))

contrasts_C13.df = data.frame(posthoc.mod3$contrasts) %>%
  tidyr::separate(contrast, into=c("tillage_1", "day_1", "tillage_2", "day_2"), remove=FALSE) %>%
  filter(tillage_1 != tillage_2, day_1 == day_2) %>%
  mutate(day = as.numeric(as.character(day_1))) %>%
  dplyr::select(-day_1, -day_2) %>%
  mutate(contrast = paste(tillage_1, tillage_2, sep="-")) %>% #
  dplyr::rename(sample_day = day) %>%
  mutate(adj.p.value = p.adjust(p.value, method = "BH"))

plot_C13_rate <- filter(contrasts_C13.df, adj.p.value<0.1) %>%
  mutate(label = ifelse(adj.p.value  <0.05, "*", ".")) %>%
  mutate(mean_min = 0.0032, Land_Management = NA, Treatment="Cellulose") %>%
  dplyr::select(mean_min, Land_Management, sample_day, label, Treatment) #sample day 6 is significant

#Xylose
df4 <- data %>%
  filter(!grepl("H2O", sample)) %>%
  filter(co2 == "45") %>%
  filter(Treatment == "13X") %>%
  mutate(Treatment = factor(Treatment),
         Rep = factor(Rep),
         sample_day = factor(sample_day)) 

mod4 <- lmer(C_rate_mg ~ Land_Management*sample_day  + (1|Rep), data=df4)
anova(mod4) #day and the interaction btw land management and day is significant (but not land management in itself)

# Run a posthoc pairwise comparisons and get the means
#lsmeans(treat.model, pairwise~ecosystem*day, adjust="tukey")
posthoc.mod4 = lsmeans(mod4, pairwise~Land_Management*sample_day, adjust=NULL)

rate_X13.df = data.frame(posthoc.mod4$lsmeans) %>%
  mutate(sample_day = as.numeric(as.character(sample_day)))

contrasts_X13.df = data.frame(posthoc.mod4$contrasts) %>%
  tidyr::separate(contrast, into=c("tillage_1", "day_1", "tillage_2", "day_2"), remove=FALSE) %>%
  filter(tillage_1 != tillage_2, day_1 == day_2) %>%
  mutate(day = as.numeric(as.character(day_1))) %>%
  dplyr::select(-day_1, -day_2) %>%
  mutate(contrast = paste(tillage_1, tillage_2, sep="-")) %>% #
  dplyr::rename(sample_day = day) %>%
  mutate(adj.p.value = p.adjust(p.value, method = "BH"))

plot_X13_rate <- filter(contrasts_X13.df, adj.p.value<0.1) %>%
  mutate(label = ifelse(adj.p.value  <0.05, "*", ".")) %>%
  mutate(mean_min = 0.022, Land_Management = NA, Treatment="Xylose") %>%
  dplyr::select(mean_min, Land_Management, sample_day, label, Treatment) #sample day 2 is significant

#### Analysis 3: Does cumulative 13C mineralization differ between tillage for each substrate/day? ####
#Perform analysis individually for each substrate
#Cellulose
mod5 <- lmer(cum_C ~ Land_Management*sample_day  + (1|Rep), data=df3)
anova(mod5) #land management, day, and their interaction is significant

# Run a posthoc pairwise comparisons and get the means
posthoc.mod5 = lsmeans(mod5, pairwise~Land_Management*sample_day, adjust=NULL)

cum_C13.df = data.frame(posthoc.mod5$lsmeans) %>%
  mutate(sample_day = as.numeric(as.character(sample_day)))

contrasts_cum_C13.df = data.frame(posthoc.mod5$contrasts) %>%
  tidyr::separate(contrast, into=c("tillage_1", "day_1", "tillage_2", "day_2"), remove=FALSE) %>%
  filter(tillage_1 != tillage_2, day_1 == day_2) %>%
  mutate(day = as.numeric(as.character(day_1))) %>%
  dplyr::select(-day_1, -day_2) %>%
  mutate(contrast = paste(tillage_1, tillage_2, sep="-")) %>% #
  dplyr::rename(sample_day = day) %>%
  mutate(adj.p.value = p.adjust(p.value, method = "BH"))

plot_cum_C13 <- filter(contrasts_cum_C13.df, adj.p.value<0.1) %>%
  mutate(label = ifelse(adj.p.value  <0.05, "*", ".")) %>%
  mutate(mean_cumC = 3, Land_Management = NA, Treatment="Cellulose") %>%
  dplyr::select(mean_cumC, Land_Management, sample_day, label, Treatment) #sample day 2 is significant

#Xylose
mod6 <- lmer(cum_C ~ Land_Management*sample_day  + (1|Rep), data=df4)
anova(mod6) #land management, day, and their interaction is significant

# Run a posthoc pairwise comparisons and get the means
posthoc.mod6 = lsmeans(mod6, pairwise~Land_Management*sample_day, adjust=NULL)

cum_X13.df = data.frame(posthoc.mod6$lsmeans) %>%
  mutate(sample_day = as.numeric(as.character(sample_day)))

contrasts_cum_X13.df = data.frame(posthoc.mod6$contrasts) %>%
  tidyr::separate(contrast, into=c("tillage_1", "day_1", "tillage_2", "day_2"), remove=FALSE) %>%
  filter(tillage_1 != tillage_2, day_1 == day_2) %>%
  mutate(day = as.numeric(as.character(day_1))) %>%
  dplyr::select(-day_1, -day_2) %>%
  mutate(contrast = paste(tillage_1, tillage_2, sep="-")) %>% #
  dplyr::rename(sample_day = day) %>%
  mutate(adj.p.value = p.adjust(p.value, method = "BH"))

plot_cum_X13 <- filter(contrasts_cum_X13.df, adj.p.value<0.1) %>%
  mutate(label = ifelse(adj.p.value  <0.05, "*", ".")) %>%
  mutate(mean_cumC = 2.1, Land_Management = NA, Treatment="Xylose") %>%
  dplyr::select(mean_cumC, Land_Management, sample_day, label, Treatment) 

#### Figure 2 ####
#plot mean C-min

c_mean <- data %>%
  filter(co2 == "45") %>%
  filter(!grepl("H2O", sample)) %>%
  group_by(co2, Treatment, sample_day, Land_Management) %>%
  summarise(mean_min = mean(C_rate_mg),
            mean_cumC = mean(mg_cum_C),
            sd_min = sd(C_rate_mg),
            sd_cumC = sd(mg_cum_C)) %>%
  ungroup() %>%
  unique() %>%
  mutate(Treatment = ifelse(Treatment == "12C", "12C Control",
                            ifelse(Treatment == "13C", "Cellulose", "Xylose"))) %>%
  filter(Treatment != "12C Control")

#cumulative C plot
cellulose <- c_mean %>% filter(Treatment == "Cellulose")

cell_cum_plot <- ggplot(cellulose, aes(x=sample_day,y=mean_cumC, color=Land_Management)) +
  geom_point() +
  theme_linedraw() +
  geom_line() +
  geom_errorbar(aes(ymin=mean_cumC-sd_cumC, ymax=mean_cumC+sd_cumC), width=0.2) +
  facet_wrap(~Treatment) +
  scale_color_npg(labels=c("No till", "Plow till"), name = "Tillage") +
  xlab("Time since substrate addition (days)") +
  labs(y=expression("Cumulative "^{13}*"C mineralization (mg)")) +
  theme(strip.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        legend.position = c(0.9, 0.15),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  geom_text(data=plot_cum_C13, aes(x=sample_day, y=mean_cumC, label=label), size=6, inherit.aes=FALSE) 
# geom_text(data=plot_cum_X13, aes(x=sample_day, y=mean_cumC, label=label), size=6, inherit.aes=FALSE)
cell_cum_plot

xylose <- c_mean %>% filter(Treatment == "Xylose")

xyl_cum_plot <- ggplot(xylose, aes(x=sample_day,y=mean_cumC, color=Land_Management)) +
  geom_point() +
  theme_linedraw() +
  geom_line() +
  geom_errorbar(aes(ymin=mean_cumC-sd_cumC, ymax=mean_cumC+sd_cumC), width=0.2) +
  facet_wrap(~Treatment) +
  scale_color_npg(labels=c("No till", "Plow till"), name = "Tillage") +
  xlab("Time since substrate addition (days)") +
  labs(y=expression("Cumulative "^{13}*"C mineralization (mg)")) +
  theme(strip.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        legend.position = "none") +
  geom_text(data=plot_cum_X13, aes(x=sample_day, y=mean_cumC, label=label), size=6, inherit.aes=FALSE)
xyl_cum_plot

plot_grid(cell_cum_plot + theme(axis.title.x = element_blank()), xyl_cum_plot, 
          nrow=2, labels = "AUTO", align = "hv")

##### Figure 3 ####
cell_min_plot <- ggplot(cellulose, aes(x=sample_day,y=mean_min, color=Land_Management)) +
  geom_point() +
  theme_linedraw() +
  geom_line() +
  geom_errorbar(aes(ymin=mean_min-sd_min, ymax=mean_min+sd_min), width=0.2) +
  facet_wrap(~Treatment, nrow=1, scales="free") +
  scale_color_npg(labels=c("No till", "Plow till"), name="Tillage") +
  xlab("Time since substrate addition (days)") +
  labs(y=expression(""^{13}*"C mineralization rate (mg/hr)")) +
  theme(strip.text.x = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        legend.position = c(0.9, 0.4),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  labs(y=expression(""^{13}*"C mineralization rate (mg/hr)")) +
  geom_text(data=plot_C13_rate, aes(x=sample_day, y=mean_min, label=label), size=6, inherit.aes=FALSE) 
cell_min_plot

xyl_min_plot <- ggplot(xylose, aes(x=sample_day, y=mean_min, color=Land_Management))+
  geom_point() +
  theme_linedraw() +
  geom_line() +
  geom_errorbar(aes(ymin=mean_min-sd_min, ymax=mean_min+sd_min), width=0.2) +
  facet_wrap(~Treatment, nrow=1, scales="free") +
  scale_color_npg(labels=c("No till", "Plow till"), name="Tillage") +
  xlab("Time since substrate addition (days)") +
  labs(y=expression(""^{13}*"C mineralization rate (mg/hr)")) +
  theme(strip.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.position = "none")+
  labs(y=expression(""^{13}*"C mineralization rate (mg/hr)")) +
  geom_text(data=plot_X13_rate, aes(x=sample_day, y=mean_min, label=label), size=6, inherit.aes=FALSE)
xyl_min_plot

plot_grid(cell_min_plot, xyl_min_plot, 
          nrow=2, labels = "AUTO", align = "hv")

