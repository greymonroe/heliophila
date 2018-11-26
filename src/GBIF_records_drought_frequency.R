#Heliophila data cleaning
#9/6/2018

#accompanies/partially replicates JGM script "GBIF_records_drought_frequency.R" rows 1-19
#things to clean/check:
#basis of specimen
#names through TNRS
#duplicates?

library(raster)
library(tidyverse)
library(ggplot2)
#library(plotly)
#library(rgbif) #0.9.8
#library(taxize)
# library(dplyr)

# Read occurance data 
GBIF<-read.csv("data/0006178-160311141623029/occurrence.csv");cat("\nrows in raw GBIF data file:", nrow(GBIF))
# GBIF<-read.csv("Heliophila_occurrence.csv");cat("\nrows in raw GBIF data file:", nrow(GBIF), file = "analysis log")
# summary(GBIF$species)
#261 rows have blanks for species, NA for speciesKey
# colnames(GBIF)

# select species with life history data 
species_life_history<-read.csv("data/Heliophila_life_history.csv")
species_life_history$species<-gsub("_", " ", species_life_history$species)
# unique(species_life_history$species)#42
species_life_history$species <- as.factor(species_life_history$species)
# levels(species_life_history$species) #42

# class(GBIF$species)
# [1] "factor"
# class(species_life_history$species)
# [1] "factor"
GBIF_JGMcleaned<-subset(GBIF, species %in% species_life_history$species);cat("\nrows in raw GBIF data after choosing species with life history JGM:", nrow(GBIF_JGMcleaned), file = "analysis log", append=T)
#6634
GBIF_JGMcleaned<-merge(GBIF_JGMcleaned, species_life_history)

####check names through tnrs, see how JGM cleaned and tnrs cleaned lists vary####
# # RESULT: original JGM cleaned list fine for this stage
# tnrsnames <- sapply(species_life_history$species, function(x) name_suggest(x)$key[1], USE.NAMES=FALSE)
# unique(tnrsnames)#42
# class(tnrsnames)#integer
# class(GBIF$speciesKey)#integer
# 
# GBIF_tnrscleaned<-subset(GBIF, speciesKey %in% tnrsnames);cat("\nrows in raw GBIF data after matching species keys with life history:", nrow(GBIF_tnrscleaned), file = "analysis log", append=T)
# #4546
# class(GBIF_tnrscleaned$speciesKey)#integer
# 
# unique(GBIF$speciesKey) #93, one NA
# unique(GBIF$species) #93, on blank
# 
# unique(GBIF_tnrscleaned$speciesKey) #31 no blank/NA
# unique(GBIF_tnrscleaned$species) #31 no blank/NA
# unique(GBIF_JGMcleaned$species) #42
# unique(GBIF_JGMcleaned$speciesKey) #42
# setdiff(unique(GBIF_JGMcleaned$species), unique(GBIF_tnrscleaned$species))
# # [1] "Heliophila scoparia"     "Heliophila variabilis"   "Heliophila suavissima"   "Heliophila pectinata"    "Heliophila subulata"    
# # [6] "Heliophila seselifolia"  "Heliophila trifurca"     "Heliophila glauca"       "Heliophila juncea"       "Heliophila polygaloides"
# # [11] "Heliophila macrosperma" 
# setdiff(unique(GBIF_tnrscleaned$species), unique(GBIF_JGMcleaned$species)) #none
# 
# setdiff(unique(GBIF_JGMcleaned$speciesKey), unique(GBIF_tnrscleaned$speciesKey))
# # [1] 5374503 5374416 5374488 5374414 5374353 5374397 5374495 5374383 5374453 5374307 5374341
# 
# ###matching by key and matching by name give different numbers of results because 
# # tnrs call above gives only one key per species and 11 of these species have multiple keys
# # GBIF may use one while tnrs gave me the other.
# # BUT ANYWAY, all the names in species_life_history look clean

# 
# tnrsGBIFnames <- sapply(GBIF$species, function(x) name_suggest(x)$key[1], USE.NAMES=FALSE)
# 
# gnrNames <- gnr_resolve(names=species_life_history$species, resolve_once=T, data_source_ids=11)
# #gbif backbone tax is datasource 11

####limit basis of specimen to perserved specimen (herbarium specimen)####
GBIF_JGMcleaned2<-subset(GBIF_JGMcleaned, basisOfRecord %in% "PRESERVED_SPECIMEN");cat("\nrows in raw GBIF data after choosing records from preserved speciemns:", nrow(GBIF_JGMcleaned2), file = "analysis log", append=T)

####remove records with has geospatial issues?####
GBIF_JGMcleaned3<-subset(GBIF_JGMcleaned2, !is.na(decimalLatitude));cat("\nrows in raw GBIF data after choosing records with decimal lat:", nrow(GBIF_JGMcleaned3), file = "analysis log", append=T)

# checkGSissues <- subset(GBIF_JGMcleaned3, hasGeospatialIssues==TRUE) #23 records
# droplevels(checkGSissues)
# summary(checkGSissues$issue)
# COORDINATE_ROUNDED;COORDINATE_REPROJECTED;COUNTRY_COORDINATE_MISMATCH 
# 1
# COUNTRY_COORDINATE_MISMATCH;GEODETIC_DATUM_ASSUMED_WGS84 
# 21 
# COUNTRY_COORDINATE_MISMATCH;ZERO_COORDINATE;GEODETIC_DATUM_ASSUMED_WGS84 
# 1
###can maybe save these (by georeferencing) if necesary, for now remove
GBIF_JGMcleaned4<-subset(GBIF_JGMcleaned3, hasGeospatialIssues==FALSE);cat("\nrows in raw GBIF data after removing records with geospatial issues:", nrow(GBIF_JGMcleaned4), file = "analysis log", append=T)
#2723

####remove duplicates####
#based on identical values of year, month, day, decimalLatitude, decimalLongitude, and species. 

# GBIF_JGMcleaned4$yr_lat_long <- paste0(GBIF_JGMcleaned4$year, "_",GBIF_JGMcleaned4$decimalLatitude,"_", GBIF_JGMcleaned4$decimalLongitude)
# length(unique(GBIF_JGMcleaned4$yr_lat_long))#1804
# GBIF_JGMcleaned4$yr_lat_lon_sp <- paste0(GBIF_JGMcleaned4$yr_lat_long,"_",GBIF_JGMcleaned4$species)
# length(unique(GBIF_JGMcleaned4$yr_lat_lon_sp))#2247
GBIF_JGMcleaned4$ymd_lat_lon_sp <- paste0(GBIF_JGMcleaned4$year,"_",GBIF_JGMcleaned4$month,"_",GBIF_JGMcleaned4$day,"_",
                                          GBIF_JGMcleaned4$decimalLatitude,"_", GBIF_JGMcleaned4$decimalLongitude,"_",
                                          GBIF_JGMcleaned4$species)
# length(unique(GBIF_JGMcleaned4$ymd_lat_lon_sp))#2385

# GBIF_JGMcleaned5 <- GBIF_JGMcleaned4[unique(GBIF_JGMcleaned4$ymd_lat_lon_sp),];cat("\nrows in raw GBIF data after eliminating records witn duplicate values of yr_lat_long_sp:", nrow(GBIF_JGMcleaned5), file = "analysis log", append=T)
GBIF_JGMcleaned5 <- GBIF_JGMcleaned4[!duplicated(GBIF_JGMcleaned4$ymd_lat_lon_sp),];cat("\nrows in raw GBIF data after eliminating records witn duplicate values of yr_lat_long_sp:", nrow(GBIF_JGMcleaned5), file = "analysis log", append=T)
#2385

# check1 <- GBIF_JGMcleaned4[duplicated(GBIF_JGMcleaned4[104:106]) | duplicated(GBIF_JGMcleaned4[104:106], fromLast=TRUE),]
# check1 <- GBIF_JGMcleaned4[duplicated(GBIF_JGMcleaned4$ymd_lat_lon_sp) | duplicated(GBIF_JGMcleaned4$ymd_lat_lon_sp, fromLast=TRUE),]
# checkdups <- setdiff(GBIF_JGMcleaned4,GBIF_JGMcleaned5)
# write.table(checkdups, "checkdups.txt", sep = "\t", quote = FALSE)

#May still miss duplicates that were seperately and nonidentically georeferenced, but this should be minor?
#See code from Emily Bellis below to subsample based on error radius from location:
# library(dismo)
# library(sp)
# library(rgdal)
# ###sampling bias (subsample occurrence records) for sorghum
# coordinates(sorg) <- ~lon+lat
# r <- raster(sorg) # create a RasterLayer with the extent of shgeo
# res(r) <- 0.01 #set resolution to 0.01 degree
# r <- extend(r, extent(r)+1) # extend the extent of the RasterLayer
# acsel <- gridSample(sorg, r, n=1) # sample; can also use this to split into training and test sets
# nrow(acsel) 


# Setup -------------------------------------------------------------------

library(raster)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)

GBIF_cleaned<-GBIF_JGMcleaned5
GBIF_cleaned<-subset(GBIF_cleaned, !is.na(decimalLatitude))

# identify and remove latitude outleirs

area <- rgdal::readOGR("Data/ne_10m_ocean/ne_10m_ocean.shp")

pdf("figures/Supp_records_all_map.pdf", width=12, height=6)
plot(area)
points(GBIF_cleaned$decimalLongitude, GBIF_cleaned$decimalLatitude, pch=3, col="red")
dev.off()

pdf("Figures/Supp_fig_lat_long_dist.pdf")
hist(GBIF_cleaned$decimalLatitude, breaks=100, col="gray", xlab="Heliophila specimen latitude")
abline(v=-20, lty=2)
hist(GBIF_cleaned$decimalLongitude, breaks=100, col="gray", xlab="Heliophila specimen longitude")
abline(v=40, lty=2)
abline(v=10, lty=2)
dev.off()

latitude_outliers<-which(GBIF_cleaned$decimalLatitude>-20)
GBIF_cleaned<-GBIF_cleaned[-latitude_outliers,] ;cat("\nrows in raw GBIF data after choosing species with reasonable lat:", nrow(GBIF_cleaned), file = "analysis log", append=T)

longitude_outliers<-which(GBIF_cleaned$decimalLongitude>40)
GBIF_cleaned<-GBIF_cleaned[-longitude_outliers,] ;cat("\nrows in raw GBIF data after choosing species with reasonable long::", nrow(GBIF_cleaned), file = "analysis log", append=T)

pdf("Figures/Supp_records_cleaned_map.pdf", width=30, height=20)
plot(area, xlim=c(10,40), ylim=c(-35, -20), col="gray90")
points(GBIF_cleaned$decimalLongitude, GBIF_cleaned$decimalLatitude, pch=3, col="black")

par(mfrow=c(6,7))
for(s in 1:length(unique(GBIF_cleaned$species))){
  species_s<-as.character(unique(GBIF_cleaned$species)[s])
  GBIF_speies<-subset(GBIF_cleaned, species==species_s)
  plot(area, xlim=c(10,40), ylim=c(-35, -20), col="gray90", main=species_s)
  points(GBIF_speies$decimalLongitude, GBIF_speies$decimalLatitude, pch=3, col="black")
}
dev.off()

# drought freq data
  # data files are named based on seasons in northern hemisphere, 
  # and must be reveresed because we are looking in southern hemisphere
Winter<-raster("Data/drought/Summer_drought_freq")
Spring<-raster("Data/drought/Fall_drought_freq")
Summer<-raster("Data/drought/Winter_drought_freq")
Fall<-raster("Data/drought/Spring_drought_freq")

# Extract drought frequency variables -------------------------------------

GBIF_cleaned$Winter<-raster::extract(Winter, cbind(GBIF_cleaned$decimalLongitude, GBIF_cleaned$decimalLatitude))
GBIF_cleaned$Spring<-raster::extract(Spring, cbind(GBIF_cleaned$decimalLongitude, GBIF_cleaned$decimalLatitude))
GBIF_cleaned$Summer<-raster::extract(Summer, cbind(GBIF_cleaned$decimalLongitude, GBIF_cleaned$decimalLatitude))
GBIF_cleaned$Fall<-raster::extract(Fall, cbind(GBIF_cleaned$decimalLongitude, GBIF_cleaned$decimalLatitude))

# select only records that have satellite detected drought (are located on pixels classified as land)
GBIF_cleaned<-subset(GBIF_cleaned, !is.na(GBIF_cleaned$Summer)); cat("\nrows in raw GBIF data after choosing species with climate data (found on land):", nrow(GBIF_cleaned), file = "analysis log", append=T)


# plot GBIF records and drought frequency
pdf("Figures/All_records_cleaned_map.pdf", width=10, height=10)
plot(area, xlim=c(15,35), ylim=c(-30, -25), col="white")
points(subset(GBIF_cleaned, life_history=="p")$decimalLongitude, subset(GBIF_cleaned, life_history=="p")$decimalLatitude,  cex=2, lwd=2, pch=3, col="orange2")
plot(area, xlim=c(15,35), ylim=c(-30, -25), col="white")
points(subset(GBIF_cleaned, life_history=="a")$decimalLongitude, subset(GBIF_cleaned, life_history=="a")$decimalLatitude, cex=2, lwd=2, pch=3, col="dodgerblue")

dev.off()


pdf("Figures/South_Africa_drought_GBIF.pdf", width=10, height=10)
col_scale=colorRampPalette(c("green4","green3","white","purple4","black"))(255)
seasons=c("Winter","Spring","Summer","Fall")

par(mfrow=c(2,2))
for(s in 1:length(seasons)){
  plot(get(seasons[s]), xlim=c(10,40), ylim=c(-35, -20), main=seasons[s], col=col_scale)
}
for(s in 1:length(seasons)){
  plot(get(seasons[s]), xlim=c(10,40), ylim=c(-35, -20), main=seasons[s], col=col_scale)
  points(GBIF_cleaned$decimalLongitude, GBIF_cleaned$decimalLatitude, pch=3)
}

col_scale=colorRampPalette(c("white","gray90","gray10","black"))(255)
seasons=c("Winter","Spring","Summer","Fall")

par(mfrow=c(2,2))
for(s in 1:length(seasons)){
  plot(get(seasons[s]), xlim=c(10,40), ylim=c(-35, -20), main=seasons[s], col=col_scale)
}
for(s in 1:length(seasons)){
  plot(get(seasons[s]), xlim=c(10,40), ylim=c(-35, -20), main=seasons[s], col=col_scale)
  points(GBIF_cleaned$decimalLongitude, GBIF_cleaned$decimalLatitude, pch=3)
}
dev.off()

# table of species occurance with climate data

species_counts<-table(GBIF_cleaned$species)
species_counts<-species_counts[species_counts!=0]
species_counts<-data.frame(species_counts)
colnames(species_counts)<-c("species","specimens")
write.csv(species_counts, "tables/species_counts.csv")

# Analyze drought frequency -----------------------------------------------
# calculate species means
species_means<-merge(aggregate(Winter~species,GBIF_cleaned, mean), aggregate(Spring~species,GBIF_cleaned, mean))
species_means<-merge(species_means, aggregate(Summer~species,GBIF_cleaned, mean))
species_means<-merge(species_means, aggregate(Fall~species,GBIF_cleaned, mean))
species_means<-merge(species_means, species_life_history)

write.csv(species_means, "Tables/Table_S_species_means.csv", row.names = F)

# plot drought frequency between annual and perennial
pdf("Figures/Fig_boxplots_drought_frequency_occurance.pdf", width=4, height=6)

for(s in 1:length(seasons)){
  p<-ggplot(species_means, aes(x=life_history, col=life_history, y=get(seasons[s]), label=species))+
    geom_jitter(data=GBIF_cleaned, width = .3, alpha=0.1, pch=3)+
    geom_boxplot(outlier.colour = "white", outlier.fill = "white", fill=NA, width=.3)+
    geom_jitter(data=species_means, width=0.1)+
    scale_color_manual(values=c("dodgerblue", "orange2"))+
    theme_classic()+
    theme(legend.position = "none")+
    scale_x_discrete(name="life history", labels=c("annual","perennial"))+
    labs(y="drought frequency")+
    ggtitle(seasons[s])+
    scale_y_continuous(limits = c(0,.81))
  
  plot(p)
}

for(s in 1:length(seasons)){
  p<-ggplot(species_means, aes(x=life_history, col=life_history, y=get(seasons[s]), label=species))+
    geom_boxplot(outlier.colour = "white", outlier.fill = "white", fill=NA, width=.3)+
    geom_jitter(data=species_means, width=0.1)+
    scale_color_manual(values=c("dodgerblue", "orange2"))+
    theme_classic()+
    theme(legend.position = "none")+
    scale_x_discrete(name="life history", labels=c("annual","perennial"))+
    labs(y="drought frequency")+
    ggtitle(seasons[s])+
    scale_y_continuous(limits = c(0.1,.55))
  
  plot(p)
}


dev.off()



# line plot of drought frequency over time  -------------------------------

life_history_drought_means<- species_means %>% 
  select(-species)  %>%
  gather(variable, value, -life_history) %>%
  group_by(life_history, variable) %>%
  summarise(mean = mean(value), se = sd(value)/sqrt(length(value))) %>%
  mutate(variable = reorder(variable, c(4,2,3,1)))

life_history_drought_species<- species_means %>% 
  select(-species)  %>%
  gather(variable, value, -life_history)

pdf("Figures/life_history_drought_means.pdf", height=6, width=12)

p<-ggplot(life_history_drought_means, aes(x=variable, col=life_history,group=life_history, y=mean))+
  theme_classic()+
  geom_point(size=2)+
  geom_line(lwd=1, lty=2)+
  #geom_boxplot(data=life_history_drought_species, aes(y=value))+
  geom_errorbar(width=0.25,lwd=1, aes(ymin=mean-se, ymax=mean+se))+
  scale_color_manual(values=c( "#56B4E9","#E69F00"), name="Life history", labels=c("annual", "perennial"))+
  labs(x="Season", y="Drought Frequency")
  
plot(p)
dev.off()


# analyze species means ---------------------------------------------------
install.packages(c("gee","phytools","geiger","ape","caper","MCMCglmm","phylolm","logistf"))

library(gee) 
library(phytools) 
library(geiger) 
library(ape) 
library(caper) 
library(MCMCglmm) 
library(phylolm)
library(logistf)

model_summaries<-data.frame()

# non-phylogenetic logistic regression firth penalized
species_means$life_history_num<-as.numeric(species_means$life_history)-1
row.names(species_means)<-gsub(" ", "_", species_means$species)
seasons<-c("Winter","Spring","Summer","Fall")
for(s in 1:length(seasons)){
  fit1 <- logistf(species_means$life_history_num~species_means[,seasons[s]], firth = TRUE)
  sum<-summary(fit1)
  sum_df<-data.frame(method="firth", season=seasons[s], coefficients=sum$coefficients, p=sum$prob)
  model_summaries<-rbind(model_summaries, sum_df)
}

fit1 <- logistf(species_means$life_history_num~species_means$Summer, firth = TRUE)
summary(fit1)

###Phylogenetic Logistic Regression

tree<-read.nexus("Data/TrimmedMCCwoPolytomies.nex")
tree<- rescale(tree, "depth", 1)

for(s in 1:length(seasons)){
  fit1 <-  phyloglm(species_means$life_history_num~species_means[,seasons[s]], data=species_means, phy=tree)
  sum<-summary(fit1)
  sum_df<-data.frame(method="phylo_log", season=seasons[s], coefficients=sum$coefficients[,1], p=sum$coefficients[,4])
  model_summaries<-rbind(model_summaries, sum_df)
}


write.csv(model_summaries, "Tables/Table_models.csv")

# look at dates between annual an perennial
# convert month-day to day of year, then adjust so that summer is in the middle of the year
GBIF_cleaned$date_calc<-as.Date(paste(
  GBIF_cleaned$month,
  GBIF_cleaned$day, sep="-"),"%m-%d")
GBIF_cleaned$days<-lubridate::yday(GBIF_cleaned$date_calc)
GBIF_cleaned$days_adjusted<-sapply(GBIF_cleaned$days, function(x)
  if(is.na(x)){ return(NA)
  }else if(x>135) {
    return (x-135)
  }else return (x+230))

#plot collectiond dates
pdf("Figures/collection_dates.pdf")
p<-ggplot(GBIF_cleaned, aes(x=days_adjusted, fill=life_history)) +
  geom_density(alpha=.5)+
  scale_fill_manual(values=c( "#56B4E9","#E69F00"), name="Life history", labels=c("annual", "perennial"))+
  theme_classic()+
  scale_x_continuous(name="Season", breaks=c(45,135,225, 315), labels=c("Winter","Spring","Summer","Fall"))+
  labs(title="Collection date\nof herbaria specimens", y="Density")+
  theme(plot.title = element_text(hjust = 0.5))
p
dev.off()

# TEST ANNUAL AND PERENNIAL DATES -----------------------------------------
annuals<-subset(GBIF_cleaned, life_history=="a")
perennials<-subset(GBIF_cleaned, life_history=="p")

ks.test(annuals$days,perennials$days)
t.test(days~life_history, GBIF_cleaned)
bartlett.test(days~life_history, GBIF_cleaned)

#write table of GBIF used for analyses - ie. after filtering
write.csv(GBIF_cleaned, "Tables/Table_S_GBIF_cleaned.csv")

