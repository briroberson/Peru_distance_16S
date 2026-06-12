# Load Packages ----
library(ggplot2)
library(vegan)
library(dplyr)
library(phyloseq)
library(qiime2R)
library(tidyverse)
library(lme4)
library(car)
library(bbmle)
library(lmtest)
library(ape)
library(pairwiseAdonis)
library(LDM)
library(indicspecies)
library(MASS)
library(ecole)
library(plyr)
library(mirlyn)
library(ANCOMBC)
library(ggrepel)
library(patchwork)


# Load Metadata ----

###### 1. Load data. Alternative: once you run this once, you can then save it as an R file
#and load it directly. the code for this is in section 2d and 3e

### 1a. Metadata and elevation
#the metadata
metadata<-readr::read_tsv("peru-dec24-metadata16S.tsv")

# make distance a factor
metadata$distance<- as.factor(metadata$distance)
metadata$type_ins<- as.factor(metadata$type_ins)
metadata$type<- as.factor(metadata$type)


#this is the elevation file. this is already loaded into the 16S data but can also be loaded and joined separately
#waypoints<- read.csv("F:\\Research\\waypoints.csv")

### 1b. Other files, loaded into a phyloseq
# Load Other Data ----

#load it into a phyloseq object
phy <- qza_to_phyloseq("PeruDec24_16S_table.qza", 
                       "PeruDec24_16S_rooted-tree.qza", 
                       "PeruDec24_16S_taxonomy.qza")

#add the metadata to the phyloseq because for some reason I couldn't load it directly
#order the samples in the metadata to match the order of samples in phyloseq
metadata<-metadata[ order(match(metadata$sampleID, colnames(phy@otu_table))), ]
metadata$distance<- as.factor(metadata$distance)
#add the metadata to phyloseq
phy@sam_data<- sample_data(metadata)
#make the row names of the sample data match the actual sample ID
row.names(phy@sam_data)<- metadata$sampleID

#check that metadata didn't chop off names using sample_variables(phy). this should be the header names, not data
sample_variables(phy)


# Initial Processing ----
####### 2. Filtering

### 2a. format the heading names
## Need to remove the heading from the kingdom names (they were k_Bacteria)
taxa_data<- as.data.frame(tax_table(phy)) #pull out taxa table
taxa_data$Kingdom<- gsub("^.{0,3}", "", taxa_data$Kingdom) #remove the first three characters from Kingdom column
tax_table(phy)<- as.matrix(taxa_data) #put it back into the phyloseq

#see how many mitochondria and chloroplasts there are
table(taxa_data$Family)
table(taxa_data$Order)
table(taxa_data$Kingdom)

### 2b. Filter out mitochondria, chloroplasts, and non bacteria from phyloseq object
filtered_phy <- phy %>%
  subset_taxa(   #the subset keeps rows only where the following operators are met
    Kingdom == "Bacteria" &  #this selects only where the kingdom is bacteria
      (Family  != "Mitochondria" | is.na(Family)) &  #this and the next line select things where the family and class are NOT mito/chlo (the ! means not)
      (Order   != "Chloroplast"| is.na(Order))
  )

phy
filtered_phy

### 2c. Filter out singletons 
pruned_filtered_phy<-prune_taxa(taxa_sums(filtered_phy)>1, filtered_phy) #this is similar to subset but it keeps only the taxa that had more than 1 occurence using the taxa sums function
pruned_filtered_phy
filtered_phy
#I printed them to compare and only 2 taxa was a singleton so it only removed one

### 2d. Filter out the pos and neg controls AND the 3 samples that we chose to drop based on # of ASVs because 3 sites were sampled twice for quality
final_filtered_phy<-subset_samples(pruned_filtered_phy, 
                                   !sampleID %in% c("Dec-24-30", 'Dec-24-31', "Dec-24-46") 
                                 &   !is.na(distance)) #this removes NAs (the pos and neg controls)
final_filtered_phy #print the two to compare
pruned_filtered_phy

#save it as R file so it can be easily loaded. at this point I recommend continuing
#through the rarefying step and save the rarefied file instead
saveRDS(final_filtered_phy, file = "final_filtered_phy.rds") #use whatever file path for where you want to save it
final_filt_phy<-readRDS("final_filtered_phy.rds")


######## 3. Rarefying. if you want to skip to the rarefying and not do all of these steps (3a-3d) 
#that visualize the number of reads and determine the rarefy value, skip to section 3e

###3a. Plot histogram of number of reads per sample
reads_per_sample<- data.frame(sum=sample_sums(final_filtered_phy))
# reads_per_sample$sampleID<- row.names(reads_per_sample)
# reads<-left_join(metadata, reads_per_sample, by='sampleID')
# ggplot(reads_per_sample, aes(x=sum))+
#   geom_histogram(binwidth=2500)
# 
# reads_inslat<- reads %>% 
#   filter(type_ins=='ins_latrine')
# reads_ouslat<- reads %>% 
#   filter(type_ins=='ous_latrine')
# reads_insveg<- reads %>% 
#   filter(type_ins=='ins_l=vegpatch')
# 
# ### 3b. Determine the minimum number of reads
# smin<- min(sample_sums(final_filtered_phy))
# smin #this was used in tutorials as the value to rarefy to but obviously 0 won't work
# 
# 
# ### 3c. Plot the rarefaction curves to determine sampling depth
# 
# #For all the data
# otu.matrix = otu_table(final_filtered_phy) #make data into data frame
# otu.matrix = as.data.frame(t(otu.matrix))
# 
# #plot
# otu.rarecurve = rarecurve(otu.matrix, step = 50, label = F)
# 
# #for just inside latrine
# inslat_phy<-subset_samples(final_filtered_phy, type_ins=='ins_latrine')
# inslat_otu.matrix = otu_table(inslat_phy) #make data into data frame
# inslat_otu.matrix = as.data.frame(t(inslat_otu.matrix))
# rarecurve(inslat_otu.matrix, step=50, label=F)
# 
# #just outside latrine
# ouslat_phy<-subset_samples(final_filtered_phy, type_ins=='out_latrine')
# ouslat_otu.matrix = otu_table(ouslat_phy) #make data into data frame
# ouslat_otu.matrix = as.data.frame(t(ouslat_otu.matrix))
# rarecurve(ouslat_otu.matrix, step=50, label=F)
# 
# #just veg patch
# veg_phy<-subset_samples(final_filtered_phy, type_ins=='ins_veg_patch')
# veg_otu.matrix = otu_table(veg_phy) #make data into data frame
# veg_otu.matrix = as.data.frame(t(veg_otu.matrix))
# rarecurve(veg_otu.matrix, step=50, label=F)
# 
# 
# abline(v=6038)
# abline(v=13334)
# abline(v=36055)
# 
# 
# 
# #we have looked at the curves and discussed with KBK about what value to use to rarefy
# # she suggested we compare the lower value 13334 with the higher value 36055 to see if losing 
# # that extra sample makes a difference and the results seem to be similar so we can move forward with 
# # the higher value
# 
# ### 3d. Do the rarefying using the chosen value. rngseed sets the seed for us within the function
# filt_rare_phy<- rarefy_even_depth(final_filt_phy, rngseed=200, sample.size=36055)
# filt_rare_phy
#compare the two phyloseqs just to see and confirm that the expected number of samples were dropped




#####rarefy using mirlyn
################do for 500 reps
#rarefy data
mirl_object_500<- mirl(final_filtered_phy, libsize=38083, set.seed=200, trimOTUs=T, replace=F, rep=500) #we chose this value with KBKs help. previously had tested lower and higher values and there was no difference so using higher

#make an empty object to put the ASV tables in
mirl_otu_500 <- vector("list", length(mirl_object_500))

#extract otu tables from each rarefied phyloseq and add to the empty object above
for (i in 1:length(mirl_object_500)){
  colnames(mirl_object_500[[i]]@otu_table) <- paste0(colnames(mirl_object_500[[i]]@otu_table))
  (mirl_otu_500[[i]] <- mirl_object_500[[i]]@otu_table)
}



#make metadata file with the correct samples (remove ones dropped during rarefying)
sample_id<- data.frame(final_filtered_phy@sam_data) 
sample_id$Samples<- row.names(sample_id)
sample_id<- sample_id %>% 
  filter(!Samples %in% c('Dec-24-44','Dec-24-22','Dec-24-23','Dec-24-21','Dec-24-24','Dec-24-6'))


sample_id <- sample_id$Samples

#make empty list for each sample
average_counts_500 <- vector("list", length(sample_id))

#give how many reps you will do
rep_500<-1:500
#make empty list to hold 5 dataframes
iter_list_500<- vector('list', length(rep_500))

#rewrite loop to select columns from each rep, then average them and put them in new otu table
for (i in 1:length(sample_id) ){
  for (j in rep_500){
    iter_list_500[[j]]<-dplyr::select(as.data.frame(mirl_otu_500[[j]]),i) #this selects each individual iteration's otu table and 
    iter_list_500[[j]]$ASVname<- row.names(iter_list_500[[j]])
  }
  
  sample_df_500<- reduce(iter_list_500[rep_500], full_join, by='ASVname')
  sample_df_500[is.na(sample_df_500)]<-0
  row.names(sample_df_500)<- sample_df_500$ASVname
  sample_df_500<- sample_df_500[,c(1, 3:(1+length(rep_500)))]
  sample_average_500 <- data.frame(rowMeans(sample_df_500))
  colnames(sample_average_500) <- sample_id[[i]]
  average_counts_500[[i]] <- sample_average_500
}
average_count_df_500 <- do.call(cbind, average_counts_500)

write.csv(x=average_count_df_500, file="D:\\Soil\\Dec24_16S\\500rep_averaged_OTUtable.csv")
write.csv(x=mirl_object_500, file="D:\\Soil\\Dec24_16S\\500rep_mirlobj.csv")

######## do some checks
#check that they all have the rarefied number of ASVs
# colSums(average_count_df_500)
# 
# #is this close to the number for the whole data frame?
# sum(iter_list_500[[1]]$`99`!=0)
# sum(average_count_df_500$`99` !=0)

#add to phyloseq
mirl_phyloseq <- final_filtered_phy
mirl_phyloseq@otu_table@.Data <- as.matrix(average_count_df_500)

rowSums(mirl_phyloseq@otu_table)==rowSums(average_count_df_500)  #should print a bunch of "TRUE"

#compare the two phyloseqs just to see and confirm that the expected number of samples are present
final_filtered_phy
mirl_phyloseq

#save to final phyloseq name used for analyses 
filt_rare_phy<- mirl_phyloseq

saveRDS(filt_rare_phy, file="D:\\Soil\\Dec24_16S\\filt_rare_phy_16s.rds") #use whatever file path for where you want to save it
filt_rare_phy<-readRDS("D:\\Soil\\Dec24_16S\\filt_rare_phy_16s.rds")




# Final phyloseq (filtered and rarefied) ----
#save the data as an R file so it doesn't have to be loaded each time.
#now when you start R, you can load the metadata and waypoints in step 1a. and skip
#steps 1b-3d
saveRDS(filt_rare_phy, file="filt_rare_phy.rds") #use whatever file path for where you want to save it
filt_rare_phy<-readRDS("filt_rare_phy_16s.rds")

#order metadata to match phyloseq
metadata<-metadata[ order(match(metadata$sampleID, colnames(filt_rare_phy@otu_table))), ]


## Alpha Diversity ----
############
############

###### 4. Diversity Analysis

## NECESSARY Calculate Diversity ----
### 4a. Calculate diversity
# 
# ### All data shannon
all_shan_div<-estimate_richness(filt_rare_phy, measures='Shannon')
#add sample ID as column for the left join
all_shan_div$sampleID<- row.names(all_shan_div)

#merge with the metadata so we can run a model and filter out the stuff already filtered out
metadata_filt<- metadata %>% 
  left_join(all_shan_div, by='sampleID') %>% 
  filter(!is.na(Shannon))

#all data richness
all_richness<- estimate_richness(filt_rare_phy, measures='Observed')
all_richness$sampleID<- row.names(all_richness)
all_richness$sampleID<- gsub(".", "-", all_richness$sampleID, fixed=TRUE)

# #merge with metadata
metadata_filt<- metadata_filt %>% 
  left_join(all_richness, by='sampleID')


## Models ----
### 4c. Run models wet data subset by soil age 

### Richness ----
#filter out veg samples
metadata_lat<- metadata_filt %>% 
  filter(type=='latrine') 

m_rich_dis<- glmer.nb(Observed~distance+(1|latrine), data=metadata_lat, na.action='na.omit')
summary(m_rich_dis)
Anova(m_rich_dis)

# test it with just inside factor
#removed random effect because of no convergence error. this is just latrines
m_rich_ins<- glm.nb(Observed~inside, data=metadata_lat, na.action='na.omit')
summary(m_rich_ins)
Anova(m_rich_ins)

#test latrine and veg patches with inside factor
m_rich_veg<- glmer.nb(Observed~type_ins+(1|latrine), data=metadata_filt, na.action='na.omit')
summary(m_rich_veg)
Anova(m_rich_veg)


### Shannon's Diversity ----
m_shan_div<-lmer(Shannon~distance+(1|latrine), data=metadata_lat, na.action='na.omit')
Anova(m_shan_div)
summary(m_shan_div)

#check model assumptions
plot(m_shan_div) #checking the linear relationship, does it look linear
plot(m_shan_div,  sqrt(abs(resid(.)))~fitted(.)) #checking variance, should be random
qqnorm(residuals(m_shan_div)) #checking normality

#compare to null model
m_shan_nullS<- lmer(Shannon~1+(1|latrine), data=metadata_filt, na.action='na.omit')
lrtest(m_shan_div, m_shan_nullS)

###### try with just inside and outside factor 
m_shan_ins<-lmer(Shannon~inside+(1|latrine), data=metadata_lat, na.action='na.omit')
Anova(m_shan_ins)
summary(m_shan_ins)

#check model assumptions
plot(m_shan_ins) #checking the linear relationship, does it look linear
plot(m_shan_ins,  sqrt(abs(resid(.)))~fitted(.)) #checking variance, should be random
qqnorm(residuals(m_shan_ins)) #checking normality

#compare to null model
m_shan_nullIns<- lmer(Shannon~1+(1|latrine), data=metadata_filt, na.action='na.omit')
lrtest(m_shan_ins, m_shan_nullIns)

### test ins, ous, and vegetation
m_shan_veg<-lmer(Shannon~type_ins+(1|latrine), data=metadata_filt, na.action='na.omit')
Anova(m_shan_veg)
summary(m_shan_veg)

#check model assumptions
plot(m_shan_veg) #checking the linear relationship, does it look linear
plot(m_shan_veg,  sqrt(abs(resid(.)))~fitted(.)) #checking variance, should be random
qqnorm(residuals(m_shan_veg)) #checking normality

# plot shannon diversity
ggplot(metadata_filt, aes(type_ins, Observed)) +
  geom_boxplot(alpha = 0.5, aes(fill=type_ins)) + #adds boxplot
  geom_jitter()+
  labs(x = NULL, y = "ASV Richness", title = "16S Richness") +
  scale_fill_manual(values=c('purple3','#74e374','cyan3'), guide='none')+ #colors the different treatments
  theme_bw() 
#  geom_jitter(aes(color=elevation))



##### 5. Beta Diversity analysis

############## This step is neccessary 

###### Beta diversity





########## BETA DIVERSITY
### NECESSARY Reroot the tree.----
# 5a. Reroot the tree
#It has to be binary but now it is not since we trimmed it
phy_tree<- phy_tree(filt_rare_phy) #put tree into an object
is.binary(phy_tree) #asking if it is binary. if false, go to next step

phy_tree(filt_rare_phy)<-multi2di(phy_tree) #fix the tree and put it back in the phyloseq
is.binary(phy_tree(filt_rare_phy)) #check if it's binary, should be true

## Permanova ----
### 5e. Permanova test----
set.seed(200) ###VERY IMPORTANT, always keep the same

filt_phy_lat<- subset_samples(filt_rare_phy, type=='latrine')
metadata_lat<- metadata_filt %>% 
  filter(type=='latrine') 

#######run permanova for distance
permanova_dis<- adonis2(distance(filt_phy_lat, method='wunifrac')~distance, data=metadata_lat, by='terms')
permanova_dis

#pairwise permanova 
permanova_pairwise(distance(filt_phy_lat, method='wunifrac'), grp=metadata_lat$distance, padj='holm')

#########distance with veg patch
metadata_filt$distance[c(25:27, 29:32)]<- 'veg'
permanova_pairwise(distance(filt_rare_phy, method='wunifrac'), grp=metadata_filt$distance, padj='holm')


#### 2 category permanova
permanova_ins<- adonis2(distance(filt_phy_lat, method='wunifrac')~inside, data=metadata_lat, by='terms')
permanova_ins

#pairwise permanova 
permanova_pairwise(distance(filt_phy_lat, method='wunifrac'), grp=metadata_lat$inside, padj='holm')

####### 3 category with vegetation
permanova_veg<- adonis2(distance(filt_rare_phy, method='wunifrac')~type_ins, data=metadata_filt, by='terms')
permanova_veg

#pairwise permanova 
permanova_pairwise(distance(filt_rare_phy, method='wunifrac'), grp=metadata_filt$type_ins, padj='holm')


## Permanova with separate 1 & 2 groups 
metadata_filt_grouped <- metadata_filt %>%
  mutate(group = case_when(
    type == "veg_patch" ~ "veg_patch",
    type == "latrine" & distance == 1 ~ "1",
    type == "latrine" & distance == 2 ~ "2",
    type == "latrine" & distance %in% 3:6 ~ "outside",
    TRUE ~ NA_character_ ))
metadata_filt_grouped$group <- as.factor(metadata_filt_grouped$group)

permanova_grouped<- adonis2(distance(filt_rare_phy, method='wunifrac')~group, data=metadata_filt_grouped, by='terms')
permanova_grouped

#pairwise permanova 
permanova_pairwise(distance(filt_rare_phy, method='wunifrac'), grp=metadata_filt_grouped$group, padj='holm')


### now for BRAY----
set.seed(200) ###VERY IMPORTANT, always keep the same

#run permanova for distance
permanova_disBray<- adonis2(distance(filt_phy_lat, method='bray')~distance, data=metadata_lat, by='terms')
permanova_disBray

#pairwise permanova 
permanova_pairwise(distance(filt_phy_lat, method='bray'), grp=metadata_lat$distance)

## 2 category permanova
permanova_insBray<- adonis2(distance(filt_phy_lat, method='bray')~inside, data=metadata_lat, by='terms')
permanova_insBray

#pairwise permanova 
permanova_pairwise(distance(filt_phy_lat, method='bray'), grp=metadata_lat$inside)

# 3 category with vegetation
permanova_vegBray<- adonis2(distance(filt_rare_phy, method='bray')~type_ins, data=metadata_filt, by='terms')
permanova_vegBray

#pairwise permanova 
permanova_pairwise(distance(filt_rare_phy, method='bray'), grp=metadata_filt$type_ins)


# PCoA Plot----
## WUNIFRAC
#get asv table and transpose for just latrine samples
asvLat<- as.data.frame(otu_table(filt_phy_lat))
tasvLat <- data.frame(t(asvLat), check.names = F)

#calculate the pcoa
pcoaLat<-cmdscale(d=distance(filt_phy_lat, method='wunifrac'), eig=T)

#retrieve species scores for it
spscorLat<-as.data.frame(wascores(x = pcoaLat$points, w = tasvLat))

#add the scores to the metadata
metadata_lat$axis01<- vegan::scores(pcoaLat)[,1]
metadata_lat$axis02<- vegan::scores(pcoaLat)[,2]

#use this function to calculate the hulls
find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hulls <- ddply(metadata_lat, "distance", find_hull)

#plot it for distance
ggplot(metadata_lat, aes(axis01, axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = distance, fill = distance), alpha = 0.1, show.legend = F) +
  geom_point(size = 3, aes(colour = distance)) +
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 18, face = "bold",color = "black"),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

#plot it for inside outside
micro.hulls <- ddply(metadata_lat, "inside", find_hull)
ggplot(metadata_lat, aes(axis01, axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = inside, fill = inside), alpha = 0.1, show.legend = F) +
  geom_point(size = 3, aes(colour = inside)) +
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 18, face = "bold",color = "black"),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))



## for veg and latrine
#get asv table and transpose for all samples
asvAll<- as.data.frame(otu_table(filt_rare_phy))
tasvAll <- data.frame(t(asvAll), check.names = F)

#calculate the pcoa
pcoaAll<-cmdscale(d=distance(filt_rare_phy, method='wunifrac'), eig=T)

#retrieve species scores for it
spscorAll<-as.data.frame(wascores(x = pcoaAll$points, w = tasvAll))

#add the scores to the metadata
metadata_filt$axis01<- vegan::scores(pcoaAll)[,1]
metadata_filt$axis02<- vegan::scores(pcoaAll)[,2]

#use this function to calculate the hulls
find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hulls <- ddply(metadata_filt, "type_ins", find_hull)

#plot it for distance
ggplot(metadata_filt, aes(axis01, axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = type_ins, fill = type_ins), alpha = 0.1, show.legend = F) +
  geom_point(size = 3, aes(colour = type_ins)) +
  scale_color_manual(labels=c('Inside Latrine','Vegetation Patch','Outside Latrine'),
                     values=c('purple1','#74e374', 'cyan2'))+
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))
# Simper ----

# run simper for Distance----
simper_dis<- simper(tasvLat, metadata_lat$distance)
simper_dis

#see the top 20
s_dis<- summary(simper_dis)
top10_1_2<-head(s_dis$`1_2`, n = 10)
top10_1_3<-head(s_dis$`1_3`, n = 10)
top10_1_4<-head(s_dis$`1_4`, n = 10)
top10_1_5<-head(s_dis$`1_5`, n = 10)

simpdis_asv1_2<- row.names(top10_1_2)
simpdis_asv1_3<- row.names(top10_1_3)
simpdis_asv1_4<- row.names(top10_1_4)
simpdis_asv1_5<- row.names(top10_1_5)

#get actual taxa info
taxa_dis <- as.data.frame(tax_table(filt_phy_lat)) #taxonomy
simper1_2_taxa<-taxa_dis[row.names(taxa_dis) %in% simpdis_asv1_2,]
simper1_3_taxa<-taxa_dis[row.names(taxa_dis) %in% simpdis_asv1_3,]
simper1_4_taxa<-taxa_dis[row.names(taxa_dis) %in% simpdis_asv1_4,]
simper1_5_taxa<-taxa_dis[row.names(taxa_dis) %in% simpdis_asv1_5,]


# run simper on Inside/Outside----
simper_ins<- simper(tasvLat, metadata_lat$inside)
simper_ins

#see the top 20
s_ins<- summary(simper_ins)
top10_ins<-head(s_ins$inside_outside, n = 10)

simpdis_asv_ins<- row.names(top10_ins)

#get actual taxa info
simper_ins_taxa<-taxa_dis[row.names(taxa_dis) %in% simpdis_asv_ins,]


# run simper on 3 categories----
simper_veg<- simper(tasvAll, metadata_filt$type_ins)
simper_veg

#see the top 20
s_veg<- summary(simper_veg)
top10_veg_insL<-head(s_veg$ins_latrine_ins_veg_patch, n = 10)
top10_veg_ousL<-head(s_veg$out_latrine_ins_veg_patch, n = 10)
top10_insL_ousL<-head(s_veg$ins_latrine_out_latrine, n = 10)

simp_asv_veg_insL<- row.names(top10_veg_insL)
simp_asv_veg_ousL<- row.names(top10_veg_ousL)
simp_asv_insL_ousL<- row.names(top10_insL_ousL)

#get actual taxa info
simper_veg_insL_taxa<-taxa_dis[row.names(taxa_dis) %in% simp_asv_veg_insL,]
simper_veg_ousL_taxa<-taxa_dis[row.names(taxa_dis) %in% simp_asv_veg_ousL,]
simper_insL_ousL_taxa<-taxa_dis[row.names(taxa_dis) %in% simp_asv_insL_ousL,]


###pcoa with simper labels 

#function to get lowest identified taxa 
bad_terms <- c(
  "Incertae_Sedis",
  "Subgroup_10")


get_lowest_tax <- function(x) {
  x <- as.character(x)
  
  for (i in rev(seq_along(x))) {
    if (!is.na(x[i]) &&
        x[i] != "" &&
        !any(sapply(bad_terms, function(b) grepl(b, x[i], ignore.case = TRUE)))) {
      return(x[i])
    }
  }
  return(NA)
}

tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
taxa_dis$label <- apply(
  taxa_dis[, tax_cols],
  1,
  get_lowest_tax)


## for veg and latrine
#get asv table and transpose for all samples
asvAll<- as.data.frame(otu_table(filt_rare_phy))
tasvAll <- data.frame(t(asvAll), check.names = F)

#calculate the pcoa
pcoaAll<-cmdscale(d=distance(filt_rare_phy, method='wunifrac'), eig=T)
#add the scores to the metadata
metadata_filt$axis01<- vegan::scores(pcoaAll)[,1]
metadata_filt$axis02<- vegan::scores(pcoaAll)[,2]

#retrieve species scores for it
spscorAll<-as.data.frame(wascores(x = pcoaAll$points, w = tasvAll))
spscorAll$ASV <- rownames(spscorAll)

spscorAll <- merge(spscorAll, taxa_dis[, "label", drop = FALSE],
                    by = "row.names", all.x = TRUE)
rownames(spscorAll) <- spscorAll$Row.names

veg_insL_df  <- spscorAll[spscorAll$ASV %in% simp_asv_veg_insL, ]
veg_ousL_df  <- spscorAll[spscorAll$ASV %in% simp_asv_veg_ousL, ]
insL_ousL_df <- spscorAll[spscorAll$ASV %in% simp_asv_insL_ousL, ]

#use this function to calculate the hulls
find_hull <- function(df) df[chull(df$axis01, df$axis02),]
micro.hulls <- ddply(metadata_filt, "type_ins", find_hull)

#plot it for distance
ggplot(metadata_filt, aes(axis01, axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = type_ins, fill = type_ins), alpha = 0.1, show.legend = F) +
  geom_point(size = 3, aes(colour = type_ins)) +
  scale_color_manual(labels=c('Inside Latrine','Vegetation Patch','Outside Latrine'),
                     values=c('purple1','#74e374', 'cyan2'))+
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  geom_text_repel(data = veg_insL_df,
                  aes(x = V1, y = V2, label = label),
                  size = 3, color = "black") +
  geom_text_repel(data = veg_ousL_df,
                  aes(x = V1, y = V2, label = label),
                  size = 3, color = "black") +
  geom_text_repel(data = insL_ousL_df,
                  aes(x = V1, y = V2, label = label),
                  size = 3, color = "black") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold",color = "black"),
    plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))




## DA----
# fixed effects must be a factor
str(filt_rare_phy@sam_data)
## trying differential abundance analysis
All_DA<-ancombc2(data = filt_rare_phy, tax_level = "Genus",
                  fix_formula = "type_ins", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.0, lib_cut = 0, s0_perc = 0.05,
                  group = "type_ins", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = F, pairwise = TRUE, dunnet = F, trend = F,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))
# inside latrine is the reference category here

#primary to view everything
res_prim = All_DA$res %>%
  mutate_if(is.numeric, function(x) round(x, 2))
saveRDS(res_prim, file='F:\\Research\\Dec24_16S\\DA_prim.rds')
res_prim <- readRDS("DA_prim.rds")

#view the structural zeros
tab_zero = All_DA$zero_ind
saveRDS(tab_zero, file='F:\\Research\\Dec24_16S\\DA_zero.rds')
tab_zero <- readRDS("DA_zero.rds")

tab_zero_insL_veg<- tab_zero %>% 
  filter((`structural_zero (type_ins = ins_latrine)`==T & `structural_zero (type_ins = ins_veg_patch)`==F)|
           `structural_zero (type_ins = ins_latrine)`==F & `structural_zero (type_ins = ins_veg_patch)`==T)

#view pairwise
res_pair<- All_DA$res_pair %>% 
  mutate_if(is.numeric, function(x) round(x, 2))
res_pair <- readRDS("DA_pair.rds")

#plot the log change in significant/passed sensitivity test for taxa between 
# Inside latrine and Vegetation Patch
res_insL_veg<- res_pair %>% 
  filter(q_type_insins_veg_patch<.05 & passed_ss_type_insins_veg_patch==T) %>% 
  dplyr::arrange(desc(lfc_type_insins_veg_patch)) %>% 
  dplyr::mutate(direct = ifelse(lfc_type_insins_veg_patch> 0, "Positive LFC", "Negative LFC"))

#make taxon and direction factors and add number index 
res_insL_veg$taxon<- factor(res_insL_veg$taxon, levels=res_insL_veg$taxon)
res_insL_veg$direct<- factor(res_insL_veg$direct, levels = c("Positive LFC", "Negative LFC"))


fig_insL_veg = res_insL_veg %>%
  ggplot(aes(x = taxon, y = lfc_type_insins_veg_patch, fill=direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_type_insins_veg_patch - se_type_insins_veg_patch, ymax = lfc_type_insins_veg_patch + se_type_insins_veg_patch), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
fig_insL_veg

#remove taxa for paper 
fig_insL_veg = res_insL_veg %>%
  ggplot(aes(x = taxon, y = lfc_type_insins_veg_patch, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black",
           position = position_dodge(width = 0.4)) +
  labs(y = "Log-fold change", x = NULL, title = "(a)") + 
  geom_errorbar(
    aes(
      ymin = lfc_type_insins_veg_patch - se_type_insins_veg_patch,
      ymax = lfc_type_insins_veg_patch + se_type_insins_veg_patch
    ),
    width = 0.2,
    position = position_dodge(0.05),
    color = "black"
  ) + 
  scale_x_discrete(labels = seq_along(unique(res_insL_veg$taxon))) +
  scale_fill_manual(values = c(
    "Positive LFC" = "purple1",
    "Negative LFC" = "#74e374"
  )) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold", color = "black"),
    plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
    legend.title = element_blank()
  )
fig_insL_veg


# Inside latrine and outside latrine
res_insL_ousL<- res_pair %>% 
  filter(q_type_insout_latrine<.05 & passed_ss_type_insout_latrine==T) %>% 
  dplyr::arrange(desc(lfc_type_insout_latrine)) %>% 
  dplyr::mutate(direct = ifelse(lfc_type_insout_latrine> 0, "Positive LFC", "Negative LFC"))

#make taxon and direction factors
res_insL_ousL$taxon<- factor(res_insL_ousL$taxon, levels=res_insL_ousL$taxon)
res_insL_ousL$direct<- factor(res_insL_ousL$direct, levels = c("Positive LFC", "Negative LFC"))

fig_insL_ousL = res_insL_ousL %>%
  ggplot(aes(x = taxon, y = lfc_type_insout_latrine, fill=direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_type_insout_latrine - se_type_insout_latrine, ymax = lfc_type_insout_latrine + se_type_insout_latrine), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank())
fig_insL_ousL


#for manuscript 

fig_insL_ousL = res_insL_ousL %>%
  ggplot(aes(x = taxon, y = lfc_type_insout_latrine, fill=direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_type_insout_latrine - se_type_insout_latrine, ymax = lfc_type_insout_latrine + se_type_insout_latrine), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", title = "(b)") + 
  scale_x_discrete(labels = seq_along(unique(res_insL_ousL$taxon))) +
  scale_fill_manual(values = c(
    "Positive LFC" = "purple1",
    "Negative LFC" = "cyan2"
  )) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(face="bold", size = 18), 
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold", color = "black"),
    plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
    legend.title = element_blank()
  )
fig_insL_ousL


#to put together 
fig_insL_veg /
  fig_insL_ousL

# Outside Latrine and Vegetation patches
res_ousL_veg<- res_pair %>% 
  filter(q_type_insout_latrine_type_insins_veg_patch<.05 & passed_ss_type_insout_latrine_type_insins_veg_patch==T) %>% 
  dplyr::arrange(desc(lfc_type_insout_latrine_type_insins_veg_patch)) %>% 
  dplyr::mutate(direct = ifelse(lfc_type_insout_latrine_type_insins_veg_patch> 0, "Positive LFC", "Negative LFC"))

#make taxon and direction factors
res_ousL_veg$taxon<- factor(res_ousL_veg$taxon, levels=res_ousL_veg$taxon)
res_ousL_veg$direct<- factor(res_ousL_veg$direct, levels = c("Positive LFC", "Negative LFC"))


fig_ousL_veg = res_ousL_veg %>%
  ggplot(aes(x = taxon, y = lfc_type_insout_latrine_type_insins_veg_patch, fill=direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_type_insout_latrine_type_insins_veg_patch - se_type_insout_latrine_type_insins_veg_patch, ymax = lfc_type_insout_latrine_type_insins_veg_patch + se_type_insout_latrine_type_insins_veg_patch), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
fig_ousL_veg

taxa_dis %>% 
  filter(Genus=='Incertae_Sedis') %>% 
  summarize()




#### libraries for heatmaps
library(plyr)
library(heatmaply)
library(gplots)
library(dendextend)
library(colorspace)
library(scales)
library(reshape)
library(funrar)
library(mctoolsr)
library(rvg)
library(officer)


### Dendrogram----
#make separate files for inside latrine in case you want to make a dendrogram for just that group
inslat_filt_phy<- subset_samples(filt_rare_phy, type_ins=='ins_latrine' )
inslat_meta<- metadata_lat %>% 
  filter(type_ins=='ins_latrine')

#dendrogram for all the data, by the 3 categories
# get the distance matrix
dist_all<- distance(filt_rare_phy, method='wunifrac')
# make a dendrogram using hclust on the distance
dend_all<- as.dendrogram(hclust(dist_all, method='ward.D2'))
#color the branches so at the first split, one split is one color and the other split is the other
dend_all<- color_branches(dend_all, k=2, col=c("cyan3","purple3"))
col_labels<- get_leaves_branches_col(dend_all)
col_labels <- col_labels[order(order.dendrogram(dend_all))]
dend_all <- set(dend_all, "labels_cex", 0.8)
dend_all<-place_labels(dend_all, as.character(metadata_filt$type_ins))
#this sets the dimensions of the plotting plane so that the dendrogram fits in it
par(mar = c(1,1,1,14))
plot_horiz.dendrogram(dend_all, side=F)


#glomerate it 
glomR <- tax_glom(filt_rare_phy, taxrank = 'Phylum')
datR <- psmelt(glomR)
datR$Phylum <- as.character(datR$Phylum)

Phylum_abundanceR <- aggregate(Abundance~Sample+Phylum, datR, FUN=sum)
Phylum_abundanceR <- cast(Phylum_abundanceR, Sample ~ Phylum)

row.names(Phylum_abundanceR)<- Phylum_abundanceR$Sample
Phylum_abundanceR<- Phylum_abundanceR[, -1]
Phylum_abundanceR<- as.data.frame(Phylum_abundanceR)


#calculate relative abundance
Phy_relabR<-make_relative(as.matrix(Phylum_abundanceR))
Phy_relabR<-data.frame(Phy_relabR, check.names=F)
#transpose it
Phy_relabRt<- data.frame((t(Phy_relabR)), check.names=F)


#plot it
all_heat<-plot_ts_heatmap(Phy_relabRt, metadata_filt, 0.01, "type_ins", colors=c('#e4daed','#61298f')) +
  theme(axis.text.y = element_text(size = 12, hjust = 0,
                                   margin = margin(c(0,-1,0,0))),
        axis.text.x = element_text(size = 10, angle = 30, hjust = 1, vjust = 1,
                                   margin = margin(c(-1.5,0,0,0))),
        plot.margin = unit(c(0.1,0.1,0.1,1), "cm")) +
  labs(title='')

#export the plot to a powerpoint to edit
fig_dml<- rvg::dml(ggobj = all_heat)

officer::read_pptx() %>%
  # add slide 
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(fig_dml, ph_location()) %>%
  # export slide 
  base::print(
    target ='D:\\Soil\\Dec24_16S\\All_Heatmap.pptx' ) #specify file path. powerpoint should not be made already, this makes it for you

