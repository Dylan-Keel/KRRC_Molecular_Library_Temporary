# Klamath molecular library Manuscript
# Cramer Fish Sciences
# Katie Karpenko
# katie.karpenko@fishsciences.net
# 10/7/2024

### Heatmap of fish detections ###

# Libraries
library(vegan)
library(tidyverse)
library(ggrepel)
library(readxl)
library(ggtext)
library(ggpubr)
library(patchwork)
library(ggh4x)
library(stringr)

## Set working directory
setwd("C:/Users/KatieKarpenko/Cramer Fish Sciences Dropbox/GIQ Projects/RESBEL-2023-Klamath eDNA/WORKING/") #update to personal!

## Import data

# import covariates and volume corrected reads (provided by Dylan Keel)
covar = read.csv("Analysis/R/data/klamath.meta.covariates.2023.csv") #update path!
flow_cor = read.csv("Analysis/R/data/all.fish.flow.cor.data.klamath.meta.2023.csv", check.names = FALSE) #update path!

## create dataframe to add levels of site in a downstream to upstream order (for plotting)
levels = data.frame("Site.Name"=unique(covar$Site.Name), 
                    "level"=c(7,22,13,44,20,21,35,4,16:18,42,43,11,12,6,8,23,
                              24,25,1,26,27:32,9,10,3,5,19,2,14,15,39:41,33,34,36:38,0),
                    "grouping"=c("Copco", "Beaver Creek", "Bogus Creek","Scott River (Control)","Fall Creek","Fall Creek","Hayden Creek","Iron Gate","Jenny Creek","Jenny Creek","Jenny Creek","Scott River (Control)","Scott River (Control)","Klamath (Below IGR)","Klamath (Below IGR)","Copco","Copco","Klamath (Below JCB)","Klamath (Below JCB)","Klamath (Below JCB)","Iron Gate","Klamath (Below JCB)","Klamath (Below JCB)","Klamath (Below JCB)","Klamath (Below JCB)","Klamath (Below JCB)","Klamath (Below JCB)","Klamath (Below JCB)","JC Boyle","JC Boyle","Iron Gate","Iron Gate","Klamath (Below CR)","Iron Gate","Scotch Creek","Scotch Creek","Scott River (Control)","Scott River (Control)","Scott River (Control)","Shovel Creek","Shovel Creek","Spencer Creek","Spencer Creek","Spencer Creek","0"),
                    "grouping_level"=c(2,10,5,15,9,9,13,1,7,7,7,15,15,4,4,2,2,11,11,11,1,11,11,11,11,11,11,11,3,3,1,1,8,1,6,6,15,15,15,12,12,14,14,14,0))
# remove Shovel3, no fish taxa detected
levels = subset(levels, Site.Name !="Shovel3")

# create vector of unique values for grouping level
grouping_level = 
  levels %>%
  select(grouping_level, grouping) %>%
  unique()

## Edit flow corrected df
# pivot longer
flow_long = 
  flow_cor %>%
  select(-1)%>%
  pivot_longer(!Site.Name, names_to="CommonName", values_to="reads")

# add option to use giga reads
flow_long$giga_reads=flow_long$reads/10^9

# create status df
# grab unique common names
CommonName = unique(flow_long$CommonName)
# set status for each name
status = c("Exotic", "Exotic", "Exotic", "Exotic", "Exotic", 
           "Native", "Native","Exotic", "Exotic", "Exotic", 
           "Exotic", "Native","Exotic", "Native", "Native", 
           "Exotic", "Native","Exotic", "Exotic", "Native", 
           "Native")
# create df of common names and status
status_df = data.frame(CommonName, status)
# check that it assigned correctly
view(status_df)

## Edit covariates df
# select 2023 only for covariates
covar2023 = subset(covar, Year == 2023)

# Merge various files for flow  corrected reads heat map
hm_df = 
  merge(flow_long, status_df, by="CommonName") %>%
  merge(covar2023, by="Site.Name") %>%
  merge(levels, by="Site.Name")

# convert 0 reads to NA
hm_df = 
  hm_df %>%
  mutate(across(c("reads", "giga_reads"), ~ na_if(.,0))) 

# set levels of site.name to orer x axis correctly. Based on levels_new order from levels table
hm_df$Site.Name = factor(hm_df$Site.Name, levels=levels$Site.Name[order(hm_df$level)], ordered=TRUE)
hm_df$grouping = factor(hm_df$grouping, levels=grouping_level$grouping[order(hm_df$grouping_level)], ordered=TRUE)

#hm_df$FigureName <- ifelse(grepl("spp.", hm_df$CommonName), "yes", "no")

#----------Heat map of reads ----------#

# labels for species names
spplist = unique(hm_df$CommonName)

# edit spp names to have scientific names display as italic in plot
spplist = c("Goldfish"="Goldfish", "Cyprinidae spp."="<i>Cyprinidae</i> spp.", 
            "Yellow Perch"="Yellow Perch", "Speckled Dace"="Speckled Dace","Golden Shiner"="Golden Shiner",   "Yellow Bullhead"="Yellow Bullhead","Largemouth Bass"="Largemouth Bass", 
            "Pumpkinseed Sunfish"="Pumpkinseed Sunfish", "Fathead Minnow"="Fathead Minnow", 
            "Brown Trout"="Brown Trout", "Bluegill Sunfish"="Bluegill Sunfish", 
            "Catostomidae spp."="<i>Catostomidae</i> spp.","Cottus spp."="<i>Cottus</i> spp.", 
            "Rainbow  Trout"="Rainbow Trout", "Coho Salmon"="Coho Salmon", "Black Crappie"="Black Crappie", 
            "Tui Chub"="Tui Chub","Ameiurus spp."="<i>Ameiurus</i> spp.", 
            "Chinook Salmon"="Chinook Salmon", "Entosphenus spp."="<i>Entosphenus</i> spp.", 
            "Green Sunfish"="Green Sunfish")

# labels for facet_grid
hm_labels_habitat=c("Reservoir"="Reservoir", "Stream"="Stream", "Native"="Native", "Exotic"="Exotic")

# plot with site names only
klamath.hm = 
  hm_df %>%
  ggplot(., aes(x=interaction(Site.Name,grouping_level), y=CommonName, fill= giga_reads)) + 
  geom_tile(colour="white", linewidth=0.25)+
  scale_y_discrete(limits=rev, labels = spplist)+
  scale_x_discrete(guide = "axis_nested")+
  scale_fill_gradient(low = "#d8e1cf", high = "#438484", na.value = "white", 
                      name="10<sup>9</sup> DNA reads/sec.")+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.title = element_markdown(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_markdown()) +
  labs(x="Location", y="Taxa Detected") +
  facet_grid(cols=vars(Habitat), rows= vars(status), scales = "free", 
             space= "free", labeller=as_labeller(hm_labels_habitat))


# remove tick marks and labels from X axis
klamath.hm.1_notic = 
  hm_df %>%
  ggplot(., aes(x=interaction(level, grouping_level), y=CommonName, fill= giga_reads)) + 
  geom_tile(colour="white", linewidth=0.25)+
  scale_y_discrete(limits=rev, labels=spplist)+
  #scale_x_discrete(guide = "axis_nested")+
  scale_fill_gradient(low = "#d8e1cf", high = "#438484", na.value = "white",
                      name="10<sup>9</sup> DNA reads/sec.", breaks=c(500,1000,1500,2000))+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "bottom",
        legend.justification = "center", 
        legend.title = element_markdown(size=10),
        legend.text = element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_markdown())+ 
  labs(y="Taxa Detected")+
  facet_grid(cols=vars(Habitat), rows= vars(status), scales = "free", 
             space= "free", labeller=as_labeller(hm_labels_habitat))

# export heatmap 
# aspect_ratio=2.5
# ggsave("Analysis/R/figs/klamath.hm.1.numeric_notic.png", 
#     plot=klamath.hm.1_notic, width=12, height=aspect_ratio*2.25, dpi=600)


#---------- Richness and Shannon Diversity ----------#

# subset necessary columns from larger dataframe
hm_df_subset = 
  hm_df %>%
  select(Site.Name,grouping_level,grouping, Habitat, level)

# native and exotic richness and diversity per site separately
shan.rich = 
  flow_long %>%                              #grab original df
  left_join(status_df, by="CommonName") %>%  #add status for each species
  filter(status=="Native") %>%
  group_by(Site.Name) %>%            #group by site name and status
  summarise(richness = specnumber(reads),
            shannon = diversity(reads, , index="shannon")) %>%# calculate diversity and richness
  mutate(shan.rounded = round(shannon, 2)) %>%
  left_join(hm_df_subset, by="Site.Name") %>% #add some additional data for plotting 
  mutate(across(c(richness,shannon,shan.rounded), ~ na_if(.,0))) %>%
  unique() %>%
  merge(data.frame(shan.name = "Native Species Diversity",
                   rich.name = "Native Species Richness"))


# set levels of site.name to order x axis correctly. Based on levels_new order from levels table
#shan.rich$grouping = factor(shan.rich$grouping, levels=grouping_level$grouping[order(shan.rich$grouping_level)], ordered=TRUE)

### heatmap of shannon diversity
shan.plot = 
  shan.rich %>%
  ggplot(., aes(x=interaction(level, grouping_level), y=shan.name, fill= shannon)) + 
  geom_tile(colour="black", linewidth=0.2)+
  #  geom_text(aes(label = shan.rounded), color ="black", size=1.75) +  
  scale_y_discrete(limits=rev)+
  #  scale_x_discrete(guide = "axis_nested")+
  scale_fill_gradient(low = "#C5FFC2", high = "#196B24", na.value = "white", 
                      name = "Native Species Diversity", limits=c(0,2), breaks=c(0.5,1,1.5,2))+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "none",
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),) + 
  facet_grid(cols=vars(Habitat),  scales = "free", 
             space= "free")


### heatmap of richness
rich.plot = 
  shan.rich %>%
  ggplot(., aes(x=interaction(level, grouping_level), y=rich.name, fill= richness)) + 
  geom_tile(colour="black", linewidth=0.2)+
  # geom_text(aes(label = richness), color="black", size=2) +  
  scale_y_discrete(limits=rev)+
  scale_x_discrete(guide = "axis_nested")+
  scale_fill_gradient(low = "#cdcbcb", high = "#696969", na.value = "white", 
                      name="Native Species Richness", limits=c(0,8), breaks=c(2,4,6,8))+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "none",
        legend.position = "bottom",
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.ticks.x=element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text( vjust = 1)) + 
  labs(x="Location") +
  facet_grid(cols=vars(Habitat),  scales = "free", 
             space= "free")

# Combine heatmap of reads, shannon diversity index and richness into one plot
hm.rich.shan =
  (klamath.hm.1_notic / plot_spacer() /shan.plot / plot_spacer()/rich.plot) +  
  plot_layout(widths = c(12,12, 12,12, 12), heights = unit(c(3.75,-0.3, .2, -0.3, .2), c('in', 'in')))  + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# Export combined plot. Remaining edits made in graphics editor.
aspect_ratio=2.5
ggsave("Analysis/R/figs/klamath.heatmap.richshan_notext3.png", 
       plot=hm.rich.shan, width=9.5, height=aspect_ratio*2.5, units="in", dpi=600)

#---------END-------------#