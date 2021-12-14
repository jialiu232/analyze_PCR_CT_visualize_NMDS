# @jialiu
# After getting the CT threshold from `get_CT_threshold_biomarkAOB.R`,
# I am visualizing Ct values with heat map and NMDS in this script

library(dplyr)
library(ggplot2)
library(tidyverse)
library(viridis)
library(vegan)
## Read in the data
setwd("/Users/liujia/Desktop/data_anal_nmds_anosim/")
f <- readxl::read_excel("biomark_results_sheet.xlsx")


## Data preprocessing
filt_f <- f %>% 
  filter(!startsWith(Sample, "16") & !startsWith(Sample, "amo") & !startsWith(Assay, "16")) %>% 
  mutate(Ct = ifelse(Ct > 23, 40, Ct)) %>% 
  #mutate(Sample_site = substr(Sample, 1, 1)) %>% 
  mutate(Temp = Sample) %>% 
  separate(Temp, c("Sample_site", "Plant_type", "Year", "P_var", "Nitrogen", "Duplicate", NA), "_") %>% 
  filter(!grepl("_", Assay)) %>% 
  arrange(Plant_type, Nitrogen, Year, Sample_site, Sample, Assay)

# need to simplify values in Sample 
# check Sample names by:
#df_ll <- as.data.frame(table(ll))
#dim(df_ll)
#df_ll$ll
# thus, there are 36 unique samples; each of them has
# 2 replicates; each replicate has 78 rows
filt_f$Uniq_sampleID <- c(rep(1:16, each=78*2), rep(1:20, each=78*2))
filt_f$Reps <- rep(1:2, each=78, times = 36)

filt_f <- filt_f %>% 
  mutate(Plant_type = ifelse(Plant_type=="C", "corn", "miscanthus")) %>% 
  mutate(Uniq_sampleID = sprintf("%02d", as.numeric(Uniq_sampleID))) %>% 
  unite("newSample", Plant_type, Nitrogen, Uniq_sampleID, Reps, sep= "_", remove = FALSE) %>% 
  mutate(Prim_order = substring(Assay, 7)) %>% 
  mutate(Prim_order = paste("amoA_AOB_p", Prim_order, sep = "")) %>% 
  filter(Plant_type == "corn")
  
# Assay will be renamed to Prim_order: amoA_AOB_p01 - amoA_AOB_p78

## Make heatmap of (Sample) newSamle vs (Assay) Prim_order
filt_f %>% 
  ggplot(aes(x=Prim_order, y=forcats::fct_rev(newSample), fill=Ct)) +
  geom_tile(colour="gray20", size=0.5, stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 5.5)) +
  #scale_fill_gradient(low = "black", high = "grey") +
  scale_fill_viridis(option="A", direction = -1, begin = 0, limits=c(15, 40), breaks=seq(15,40,by=5)) +
  #scale_fill_continuous() +
  xlab("Primer set") + 
  ylab("Soil sample")
  






## NMDS to PCR results
# To make the NMDS by `metaMDS, we need to first reshape "Sample", "Assay", and "Ct" into a matrix
library(reshape2)
filt_f_matrix <- filt_f %>% 
  acast(Sample~Prim_order, value.var = "Ct")
set.seed(100)

# Run NMDS
pcr_NMDS <- metaMDS(filt_f_matrix, distance = "bray", k=2, trymax = 1000)

pcr.spp.fit <- envfit(pcr_NMDS, filt_f_matrix, permutations = 999) # this fits species vectors

# Show stress plot
stressplot(pcr_NMDS)

# Plot pcr_NMDS

# 1. extract the scores (the x and y coordinates of the site (rows) 
#     and species (cols) ) and add the variables we wanted

# make the site score dataframe
data.scores <- as.data.frame(scores(pcr_NMDS, display = "sites"))
data.scores$site <- rownames(data.scores)

data.scores <- data.scores %>% 
  separate(site, c("S_site", "Plant_type", "Year", "P_var", "Nitrogen", "Duplicate", NA), "_") %>% 
  mutate(site = substr(S_site, 1, 1)) %>% 
  mutate(Plant_type = ifelse(Plant_type=="C", "Corn", "Miscanthus"))
  #mutate(Nitrogen = as.numeric(substring(Nitrogen, first = 2)))
  #mutate(Year = as.integer(Year))

head(data.scores)

# make the species score dataframe
species.scores <- as.data.frame(scores(pcr_NMDS, "species"))
species.scores$species <- rownames(species.scores)
head(species.scores)

# use env to add arrows into figure for significant primers
spp.scrs <- as.data.frame(scores(pcr.spp.fit, display = "vectors"))
spp.scrs$Species <- rownames(spp.scrs)
spp.scrs$pval <- pcr.spp.fit$vectors$pvals
spp.scrs$r2 <- pcr.spp.fit$vectors$r
sig.spp.scrs <- spp.scrs %>%  # choose the species vectors with the smallest p-values and r2 bigger than or equal to 0.6
  filter(pval <= 0.05) %>% 
  arrange(pval, desc(r2)) %>% 
  filter(r2 >= 0.6)
sig.spp.vec <- sig.spp.scrs$Species
sig.species.scores <- species.scores %>% 
  filter(species %in% sig.spp.vec)
  
# plot NMDS
ggplot() +
  geom_text(data = species.scores, aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, shape=Plant_type, color=Nitrogen), size = 6) +
  ggtitle("NMDS for Ct values by 4 nitrogen levels") +
  scale_size_manual(values=c(1.5,2.5,4,6)) +
  scale_color_viridis(option="cividis", discrete = TRUE) +
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
    axis.text.y = element_blank(), # remove y-axis text
    axis.ticks = element_blank(),  # remove axis ticks
    axis.title.x = element_text(size=18), # remove x-axis labels
    axis.title.y = element_text(size=18), # remove y-axis labels
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),  #remove major-grid labels
    panel.grid.minor = element_blank(),  #remove minor-grid labels
    plot.background = element_blank()
  ) +
  annotate(geom="text", x=-0.1, y=0.13, label="Stress = 0.0954893",
           color="black") +
  geom_segment(data = sig.species.scores, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.species.scores, aes(x=NMDS1, y=NMDS2, label = species), cex = 3, direction = "both", segment.size = 0.25) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap


## Another way to plot NMDS with a hull 
# plot pcr_NMDS in another way: plot a hull around each of the nitrogen groups
hull.data <- data.scores %>% 
  group_by(Plant_type) %>% 
  slice(chull(NMDS1, NMDS2))

ggplot() +
  geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS2, fill=factor(Plant_type)), alpha=0.5) +
  geom_text(data = species.scores, aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +
  geom_point(data = species.scores, aes(x=NMDS1,y=NMDS2),alpha=0.5) +
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, size = Year, shape=Plant_type, color=Nitrogen)) +
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Year),size=2,vjust=0) +
  #scale_color_gradient(low = "orange", high = "purple") + 
  #scale_color_gradientn(colours = rainbow(5)) 
  scale_size_manual(values=c(1.5,2.5,4,6)) +
  scale_color_viridis(option="cividis", discrete = TRUE) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank()
  )


# Some patterns can be see in this figure, but not quite clear. 
# Next we will plot under specific conditions, to see if any pattern can be see

## Within one plant type, plot Nitrogen and Year
# Within Corn, plot nitrogen
corn.data.scores <- data.scores %>% 
  filter(Plant_type=="Corn")

hull.corn.data <- corn.data.scores %>% 
  group_by(Nitrogen) %>%
  slice(chull(NMDS1, NMDS2))

ggplot() +
  #stat_ellipse(data=corn.data.scores,aes(x=NMDS1,y=NMDS2, color=Nitrogen)) +
  geom_polygon(data = hull.corn.data, aes(x=NMDS1, y=NMDS2, fill=Nitrogen), alpha=0.4) +
  #geom_text(data = species.scores, aes(x=NMDS1,y=NMDS2,label=species),alpha=0.7, size = 3) +
  geom_point(data=corn.data.scores,aes(x=NMDS1,y=NMDS2, color=Nitrogen), size = 5) +
  #scale_size_manual(values=c(1.5,2.5,4,6)) +
  scale_color_viridis(option="cividis", discrete = TRUE) +
  scale_fill_viridis_d(option="viridis", begin = 0,
                       end = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank()
  ) +
  annotate(geom="text", x=-0.12, y=0.13, label="Stress = 0.0954893",
           color="black") +
  coord_cartesian(xlim = c(-0.16, 0.19), ylim = c(-0.10, 0.14)) +
  geom_point(data=sig.species.scores,aes(x=NMDS1,y=NMDS2), size = 1) +
  geom_segment(data = sig.species.scores, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.15, "cm")), colour = "grey10", lwd=0.2) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.species.scores, aes(x=NMDS1, y=NMDS2, label = species), cex = 3, direction = c("both", "y", "x"), segment.size = 0.6) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap



# Save ggplot figure
# 1. Method 1
#png("nmds_corn_nitrogen.png")
#print(nmds_corn)
#dev.off()

# 2. Method 2. use ggsave
#ggsave("nmds_corn_2.png", plot = nmds_corn)

## ANOSIM test for corn data

# 1. filter out the similarity matrix of corn to be used for ANOSIM
filt_corn_matrix <- filt_f_matrix[which(sapply(strsplit(rownames(filt_f_matrix), "_"), `[`, 2) == "C"), ]

# 2. Get the vector "nitro" of Nitrogen level corresponding to each corn sample in filt_corn_matrix
site_N <- data.frame(Plant_type = filt_f$Plant_type, Sample = filt_f$Sample, Nitrogen = filt_f$Nitrogen)
site_N <- site_N %>%  # site_N is a temporary dataframe stores corn Sample name and their corresponding Nitrogen level
  filter(Plant_type == "corn") %>% 
  unique()

# getN is a helper method to get the corresponding Nitrogen level for given Sample name from site_N dataframe
getN <- function(x) {
  nitrogen <- site_N$Nitrogen[site_N$Sample == x]
  return(nitrogen)
}

nitro <- as.vector(sapply(rownames(filt_corn_matrix), getN))

#sampName <- sapply(strsplit(rownames(filt_corn_matrix), "_"), `[`, 1)

# 3. run ANOSIM test for corn at different nitrogen level
filt_corn_matrix.dist <- vegdist(filt_corn_matrix)  # Make dissimilarity matrix by vegdist
filt_corn_matrix.ano <- anosim(filt_corn_matrix.dist, nitro, distance = "bray", permutations = 9999)

summary(filt_corn_matrix.ano)
# plot ano
par(mar = c(4.5, 4.5, 2, 2)) # Set the margin on all sides to 2
plot(filt_corn_matrix.ano,
     xlab="Classes", ylab="Dissimilarity ranks")









# Within Miscanthus, plot nitrogen and year; you need to change between "color=Year" and "color=Nitrogen"
miscanthus.data.scores <- data.scores %>% 
  filter(Plant_type=="Miscanthus")

hull.data2 <- miscanthus.data.scores %>% 
  group_by(Nitrogen) %>% 
  slice(chull(NMDS1, NMDS2))

ggplot() +
  geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS2, fill=factor(Plant_type)), alpha=0.5) +
  geom_text(data = species.scores, aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +
  geom_point(data = species.scores, aes(x=NMDS1,y=NMDS2),alpha=0.5) +
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, size = Year, shape=Plant_type, color=Nitrogen))
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Year),size=2,vjust=0) +

# nitrogen
ggplot() +
  geom_text(data = species.scores, aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +
  geom_point(data=miscanthus.data.scores,aes(x=NMDS1,y=NMDS2, color=Nitrogen), size = 5) +
  ggtitle("NMDS for Miscanthus Ct values by Sample and Assay at different nitrogen level") +
  scale_size_manual(values=c(1.5,2.5,4,6)) +
  scale_color_viridis(option="cividis", discrete = TRUE) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank()
  ) +
  annotate(geom="text", x=-0.1, y=0.13, label="Stress = 0.09315586",
           color="black")


  
  
# year
ggplot() +
  #geom_polygon(data = hull.data2, aes(x=NMDS1, y=NMDS2, fill=factor(Nitrogen)), alpha=0.5) +
  geom_text(data = species.scores, aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +
  #geom_point(data = species.scores, aes(x=NMDS1,y=NMDS2),alpha=0.3) +
  geom_point(data=miscanthus.data.scores,aes(x=NMDS1,y=NMDS2, color=Year), size = 5) +
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Year),size=2,vjust=0) +
  #scale_color_gradient(low = "orange", high = "purple") + 
  #scale_color_gradientn(colours = rainbow(5)) 
  ggtitle("NMDS for Miscanthus Ct values by Sample and Assay at different year") +
  scale_size_manual(values=c(1.5,2.5,4,6)) +
  scale_color_viridis(option="cividis", discrete = TRUE) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank()
  ) +
  annotate(geom="text", x=-0.1, y=0.13, label="Stress = 0.09315586",
           color="black")




## At nitrogen level 0 and 200, there are data for both corn and miscanthus
n0_200.data.scores <- data.scores %>% 
  filter(Nitrogen=="N0" | Nitrogen=="N200")

ggplot() +
  geom_text(data = species.scores, aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +
  geom_point(data=n0_200.data.scores,aes(x=NMDS1,y=NMDS2, color=Nitrogen, shape = Plant_type), size = 6) +
  ggtitle("NMDS for Ct values by Sample and Assay at Nitrogen N0 and N200 level") +
  scale_size_manual(values=c(1.5,2.5,4,6)) +
  scale_color_viridis(option="cividis", discrete = TRUE) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank()
  ) +
  annotate(geom="text", x=-0.1, y=0.13, label="Stress = 0.09315586",
           color="black")

