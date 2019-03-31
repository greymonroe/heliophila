## ----setup, include = FALSE----------------------------------------------
library("papaja")
library(raster)
library(tidyverse)
library(ggplot2)
library(logistf)
library(geiger) 
library(ape) 
library(phylolm)
library(gridExtra)
library(cowplot)
library(rasterVis)
library(ggtree)
library(grid)
library(png)
library(ggimage)
setwd("~/github/heliophila/manuscript/")

## ----raw_GBIF------------------------------------------------------------
GBIF<-read.csv("../data/0006178-160311141623029/occurrence.csv")

## ----Droughtdata---------------------------------------------------------
Winter<-raster("../data/drought/Summer_drought_freq")
Spring<-raster("../data/drought/Fall_drought_freq")
Summer<-raster("../data/drought/Winter_drought_freq")
Fall<-raster("../data/drought/Spring_drought_freq")

# Winter<-raster("~/VHI/droughtfreq/summer_mild")
# Spring<-raster("~/VHI/droughtfreq/fall_mild")
# Summer<-raster("~/VHI/droughtfreq/winter_mild")
# Fall<-raster("~/VHI/droughtfreq/spring_mild")

## ----cleanup-------------------------------------------------------------
species_life_history<-read.csv("../data/Heliophila_life_history.csv")
species_life_history$species<-gsub("_", " ", species_life_history$species)
species_life_history$species <- as.factor(species_life_history$species)

GBIF <- GBIF %>%
  mutate(Winter=raster::extract(Winter, cbind(decimalLongitude, decimalLatitude))) %>%
  mutate(Spring=raster::extract(Spring, cbind(decimalLongitude, decimalLatitude))) %>%
  mutate(Summer=raster::extract(Summer, cbind(decimalLongitude, decimalLatitude))) %>%
  mutate(Fall=raster::extract(Fall, cbind(decimalLongitude, decimalLatitude))) %>%
  mutate(ymd_lat_lon_sp = paste0(year,"_",month,"_",day,"_", decimalLatitude,"_", decimalLongitude,"_", species))
  
GBIF_cleaned <- GBIF %>%
  filter(species %in% species_life_history$species) %>% 
  filter(!is.na(decimalLatitude)) %>%
  filter(hasGeospatialIssues==FALSE) %>%
  filter(decimalLatitude<(-20)) %>%
  filter(decimalLongitude<40) %>%
  filter(!duplicated(ymd_lat_lon_sp)) %>%
  merge(species_life_history) %>%
  filter(!is.na(Summer))

## ----speciesmeans--------------------------------------------------------
species_means <- GBIF_cleaned %>%
  dplyr::select(species, life_history, Winter, Summer, Spring, Fall) %>%
  group_by(species, life_history) %>%
  summarise(n=n(), Winter_mean=mean(Winter), Spring_mean=mean(Spring), Summer_mean=mean(Summer), Fall_mean=mean(Fall)) %>% 
  as.data.frame()


## ----phylomodels, results="hide"-----------------------------------------
# phylogentic constraint
tree<-read.nexus("../data/TrimmedMCCwoPolytomies.nex")
tree<- rescale(tree, "depth", 1)
tree<-keep.tip(tree, tip = gsub(" ","_", as.character(species_means$species)))

species_means$life_history_num<-as.numeric(species_means$life_history)-1
row.names(species_means)<-gsub(" ", "_", species_means$species)
seasons<-c("Winter_mean", "Spring_mean", "Summer_mean","Fall_mean")

phylomodelstable<-data.frame()

for(s in 1:length(seasons)){
  fit1 <-  phyloglm(species_means$life_history_num~species_means[,seasons[s]], data=species_means, phy=tree)
  sum<-summary(fit1)
  row.names(sum$coefficients)<-NULL
  sum_df<-data.frame(Predictor=c("Intercept", gsub("_mean", " drought freq.", seasons[s])), sum$coefficients)
  phylomodelstable<-rbind(phylomodelstable, sum_df)
}

phylomodelstable <- phylomodelstable[,c( "Predictor","Estimate","p.value")]


## ----makephylogeny, eval=T, message=F, warning="hide", results="hide"----

# source("../src/ancestral state.R") #uncomment to rerun the ancestral state reconstruction. Takes a few minutes.


pdf("../figures/phylogeny.pdf", height=5.5, width=4.5)

tree <- read.nexus("../data/TrimmedMCCwoPolytomies.nex")
ladderize(tree) -> tree
tree$tip.label<-gsub("_"," ", tree$tip.label)
rescale(tree, "depth", 1) -> tree

data <- read.csv("../data/species_means_table.csv", header=T)
data$A_P_num<-as.numeric(data$LH)-1

pies<-read.csv("../data/ancestral_state.csv")
treeplot<-ggtree(tree, layout = "rectangular")+
  xlim(c(0,2))+
  geom_tippoint(size=3,col=ifelse(data$A_P_num[match(tree$tip.label,data$Species)]=="0", "orange3", "dodgerblue3"))+
  geom_tiplab(offset = .05 , size=3)+
  theme(plot.background = element_rect(fill=alpha("white",0)))

treeplot<-inset(treeplot, nodepie(pies, cols = c("X0", "X1" ), color  = c("orange3", "dodgerblue3" )), width=2.7, height=2.7)


minima <- readPNG("../pictures/minima.png")
minima <- rasterGrob(minima, interpolate=TRUE)
#https://www.gbif.org/occurrence/1099023487

deserticola <- readPNG("../pictures/deserticola.png")
deserticola <- rasterGrob(deserticola, interpolate=TRUE)
#https://www.gbif.org/occurrence/1057389408

coronopifolia<-readPNG("../pictures/coronopifolia.png")
coronopifolia <- rasterGrob(coronopifolia, interpolate=TRUE)
#https://www.gbif.org/occurrence/1099023562

ephemera<-readPNG("../pictures/ephemera.png")
ephemera <- rasterGrob(ephemera, interpolate=TRUE)
#https://www.gbif.org/occurrence/1099023490

pics<-plot_grid(minima, deserticola, coronopifolia, ephemera,
                labels = c("b","c","d","e"), ncol = 1, hjust = -1)

plot_grid(treeplot, pics, ncol=2,labels=c("a",""), rel_widths = c(1.9,0.5))


dev.off()


## ----phylogeny, fig.cap='(ref:phylogeny)', out.width = "\\textwidth", fig.pos = "!h"----
knitr::include_graphics("../figures/phylogeny.pdf", dpi = 108)

## ----make_maps_boxplots, eval=T, message=F, results="hide"---------------


pdf("../figures/maps_boxplots.pdf", width=4.5, height=7.2)
world <- map_data("world") # we already did this, but we can do it again

perennialmap<-ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group=group),  fill = NA, color = "black") + 
  coord_fixed(xlim=c(13,35), ylim=c(-35, -22)) + 
  geom_point(GBIF_cleaned %>% filter(life_history=="p"), mapping=aes(x=decimalLongitude, y=decimalLatitude), pch=3, col="dodgerblue3", alpha=0.3) +
  theme_classic()+
  theme( text = element_text(size=8))+
  labs(title="Perennial species", x="Longitude (°)", y="Latitude (°)")+
  annotate("text", x=32, y=-34, label=paste("n obs. =", nrow(GBIF_cleaned %>% filter(life_history=="p"))), size=2)

annualmap<-ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group=group),  lwd=.5, fill = NA, color = "black") + 
  coord_fixed(xlim=c(13,35), ylim=c(-35, -22)) + 
  geom_point(GBIF_cleaned %>% filter(life_history=="a"), mapping=aes(x=decimalLongitude, y=decimalLatitude), pch=3, col="orange3", alpha=0.3) +
  theme_classic()+
  theme(text = element_text(size=8))+
  labs(title="Annual species", x="Longitude (°)", y="Latitude (°)")+
  annotate("text", x=32, y=-34, label=paste("n obs. =", nrow(GBIF_cleaned %>% filter(life_history=="a"))), size=2)



col_scale<-c( 'white',"gray90", "gray10", 'black')
#col_scale<-c( "green4","green", "white","purple4","black")


Winter_SA <-crop(Winter, rbind(c(12,37),c(-35, -18)))

wintermap<-gplot(Winter_SA) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradientn(colors=col_scale, na.value="white", name=NULL) +
  geom_polygon(data = world, aes(x=long, y = lat, group=group),  fill = NA, color = "black") + 
  coord_fixed(xlim=c(13,35), ylim=c(-35, -22)) + 
  theme_classic()+
  theme(legend.position = c(0.95,.3),legend.text=element_text(size=rel(0.7)), legend.key.size=unit(.4,"line"), legend.background = element_rect(fill=alpha('white', 0)), text = element_text(size=8))+
  labs(title="Winter drought freq.", x="Longitude (°)", y="Latitude (°)")

Spring_SA <-crop(Spring, rbind(c(12,37),c(-35, -18)))

springmap<-gplot(Spring_SA) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradientn(colors=col_scale, na.value="white", name=NULL) +
  geom_polygon(data = world, aes(x=long, y = lat, group=group),  fill = NA, color = "black") + 
  coord_fixed(xlim=c(13,35), ylim=c(-35, -22)) + 
  theme_classic()+
  theme(legend.position = c(0.95,.3),legend.text=element_text(size=rel(0.7)), legend.key.size=unit(.4,"line"), legend.background = element_rect(fill=alpha('white', 0)), text = element_text(size=8))+
  labs(title="Spring drought freq.", x="Longitude (°)", y="Latitude (°)")

Summer_SA <-crop(Summer, rbind(c(12,37),c(-35, -18)))

summermap<-gplot(Summer_SA) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradientn(colors=col_scale, na.value="white", name=NULL) +
  geom_polygon(data = world, aes(x=long, y = lat, group=group),  fill = NA, color = "black") + 
  coord_fixed(xlim=c(13,35), ylim=c(-35, -22)) + 
  theme_classic()+
  theme(legend.position = c(0.95,.3),legend.text=element_text(size=rel(0.7)), legend.key.size=unit(.4,"line"), legend.background = element_rect(fill=alpha('white', 0)), text = element_text(size=8))+
  labs(title="Summer drought freq.", x="Longitude (°)", y="Latitude (°)")

Fall_SA <-crop(Fall, rbind(c(12,37),c(-35, -18)))

fallmap<-gplot(Fall_SA) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradientn(colors=col_scale, na.value="white", name=NULL) +
  geom_polygon(data = world, aes(x=long, y = lat, group=group),  fill = NA, color = "black") + 
  coord_fixed(xlim=c(13,35), ylim=c(-35, -22)) + 
  theme_classic()+
  theme(legend.position = c(0.95,.3),legend.text=element_text(size=rel(0.7)), legend.key.size=unit(.4,"line"), legend.background = element_rect(fill=alpha('white', 0)), text = element_text(size=8))+
  labs(title="Fall drought freq.", x="Longitude (°)", y="Latitude (°)")

GBIF_cleaned_long<- GBIF_cleaned %>%
  dplyr::select(species, life_history, Winter, Spring, Summer,Fall) %>%
  gather(season, drought_freq, Winter, Spring, Summer,Fall)

GBIF_cleaned_long$season<-factor(GBIF_cleaned_long$season, levels=c("Winter","Spring","Summer","Fall"))

boxplots<-ggplot(GBIF_cleaned_long, aes(x=life_history, col=life_history, y=drought_freq))+
  facet_wrap(~season, ncol=4 )+
  geom_jitter(width = .3, alpha=0.05, pch=3)+
  geom_boxplot(outlier.colour = "white", outlier.fill = NA, fill=NA, width=.3)+
  scale_color_manual(values=c("orange3","dodgerblue3"))+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=8))+
  scale_x_discrete(name="Life history", labels=c("Annual","Perennial"))+
  labs(y="Drought frequency")+
  annotate("text", x=1.5, y=0.8, label= "**")



maps<-plot_grid(annualmap, 
          perennialmap, 
          wintermap, 
          springmap,
          summermap,
          fallmap, 
          labels = "auto", ncol = 2)
plot_grid(maps, boxplots, ncol=1,labels=c("","g"), rel_heights = c(.75,.35))


dev.off()


## ----mapsboxplots, fig.cap='(ref:mapsboxplots)', out.width = "\\textwidth", fig.pos = "!h"----
knitr::include_graphics("../figures/maps_boxplots.pdf", dpi = 108)

## ----ttests, eval=T, message=F, warning="hide", results="hide"-----------
twin<-t.test(Winter~life_history, GBIF_cleaned)
tspr<-t.test(Spring~life_history, GBIF_cleaned)
tsum<-t.test(Summer~life_history, GBIF_cleaned)
tfal<-t.test(Fall~life_history, GBIF_cleaned)

melted<-GBIF_cleaned %>% select(life_history, Winter, Spring, Summer,Fall, species) %>% gather(Winter, Spring, Summer,Fall, key="season",value="droughtfreq")

library(lmerTest)
library(lme4)
fit<-lmer(droughtfreq~life_history*season+(1|species), melted)
anovatable<-data.frame(predictor=c("life history","season","life history x season"), anova(fit))

library(emmeans)
post_hoc<-emmeans(fit, list(pairwise ~ life_history*season), adjust = "tukey")
post_hoc<-as.data.frame(post_hoc$`pairwise differences of life_history, season`)
post_hoc[c(1,14,23,28), ]



## ----modelstable, results = 'asis', echo = F-----------------------------

modelstable<-bind_cols(phylomodelstable)

modelstable<- phylomodelstable %>%
  dplyr::select(Predictor, Estimate, p.value) 

colnames(modelstable)<-c("Predictor", "Estimate","P")

apa_table(modelstable
  , caption = "Phylogenetic logistic regressions between life history, and the mean drought frequency observed at specimen sites of Heliophila species the winter, spring, summer, and fall."
  , note = "Annual species were scored as 0 and perennial species as 1.",
  midrules=c(2,4,6,8),
  digits=4
)


## ----makelineplots, eval=T, message=F, results="hide", warning=F---------

pdf("../figures/line_and_dates.pdf", height=6, width=5)

life_history_drought_means<- species_means %>% 
  dplyr::select(life_history, Winter_mean, Spring_mean, Summer_mean,Fall_mean)  %>%
  gather(variable, value, -life_history) %>%
  group_by(life_history, variable) %>%
  summarise(mean = mean(value), se = sd(value)/sqrt(length(value))) %>%
  mutate(variable = reorder(variable, c(4,2,3,1)))

lines<-ggplot(life_history_drought_means, aes(x=life_history, col=life_history,group=life_history, y=mean))+
  theme_bw()+
  geom_point(size=2)+
  #geom_line(lwd=0.5, lty=1)+
  facet_grid(~variable, labeller = as_labeller(c("Winter_mean"="Winter","Spring_mean"="Spring","Summer_mean"="Summer","Fall_mean"="Fall")))+
  geom_errorbar(width=0.25,lwd=1, aes(ymin=mean-se, ymax=mean+se))+
  scale_y_continuous(limits=c(0.3,0.45), labels = c("0.300", "0.350", "0.400", "0.450"))+
  scale_x_discrete(labels=NULL)+
  scale_color_manual(values=c( "orange3","dodgerblue3"), name="Life history", labels=c("Annual", "Prennial"))+
  labs(x=NULL,y="Drought Frequency")+
  theme(axis.ticks.x =element_blank(), strip.background =element_rect(fill="white"))

GBIF_cleaned$date_calc<-as.Date(paste(
  GBIF_cleaned$month,
  GBIF_cleaned$day, sep="-"),"%m-%d")
GBIF_cleaned$days<-lubridate::yday(GBIF_cleaned$date_calc)
GBIF_cleaned$days_adjusted<-sapply(GBIF_cleaned$days, function(x)
  if(is.na(x)){ return(NA)
  }else if(x>135) {
    return (x-135)
  }else return (x+230))


dates<-ggplot(GBIF_cleaned %>% filter(is.finite(days_adjusted)), aes(x=days_adjusted, fill=life_history)) +
  geom_density(alpha=.5)+
  scale_fill_manual(values=c("orange3","dodgerblue3"), name="Life history", labels=c("Annual", "Perennial"))+
  theme_bw()+
  scale_x_continuous(name="Day of photoperiodic year", breaks=c(45,135,225, 315), labels=c("45\n(Winter)","135\n(Spring)","225\n(Summer)","315\n(Fall)"))+
  labs(title="Timing of observation for occurrence records", y="Density")+
  theme( plot.title = element_text(hjust = 0.5, size=10))

tree<-read.nexus("../data/TrimmedMCCwoPolytomies.nex")
tree<- rescale(tree, "depth", 1)
tree<-keep.tip(tree, tip = gsub(" ","_", as.character(species_means$species)))

species_means$life_history_num<-as.numeric(species_means$life_history)-1
row.names(species_means)<-gsub(" ", "_", species_means$species)

fit <-  phyloglm(species_means$life_history_num~species_means$Spring_mean, data=species_means, phy=tree)
cc <- coef(fit)
spcurve<-plogis(cc[1]+cc[2]*sort(species_means$Spring_mean))

fit <-  phyloglm(species_means$life_history_num~species_means$Summer_mean, data=species_means, phy=tree)
cc <- coef(fit)
scurve<-plogis(cc[1]+cc[2]*sort(species_means$Summer_mean))

fit <-  phyloglm(species_means$life_history_num~species_means$Fall_mean, data=species_means, phy=tree)
cc <- coef(fit)
fcurve<-plogis(cc[1]+cc[2]*sort(species_means$Fall_mean))

fit <-  phyloglm(species_means$life_history_num~species_means$Winter_mean, data=species_means, phy=tree)
cc <- coef(fit)
wicurve<-plogis(cc[1]+cc[2]*sort(species_means$Winter_mean))

curvedataframe<-data.frame(
  droughtfreq=c(sort(species_means$Winter_mean),
                sort(species_means$Spring_mean),
                sort(species_means$Summer_mean),
                sort(species_means$Fall_mean)),
  season=rep(c("Winter","Spring","Summer","Fall"), each=42),
  curves=c(wicurve, spcurve, scurve, fcurve),
  lifehistory=c(species_means$life_history_num[order(species_means$Winter_mean)],
                species_means$life_history_num[order(species_means$Spring_mean)],
                species_means$life_history_num[order(species_means$Summer_mean)],
                species_means$life_history_num[order(species_means$Fall_mean)])
                           )
curvedataframe$lifehistory<-curvedataframe$lifehistory
curvedataframe$season<-factor(curvedataframe$season, levels=c("Winter","Spring","Summer","Fall"))

fit_stats<-data.frame(season=factor(c("Winter","Spring","Summer","Fall"),levels=c("Winter","Spring","Summer","Fall")), label=c("","*","**","**"), lifehistory="1")

fitplots<-ggplot(curvedataframe, aes(x=droughtfreq, col=as.factor(lifehistory), y=lifehistory))+
  geom_jitter(height = .03, width=0, alpha=0.7)+
  facet_grid(~season)+
  theme_bw()+
  scale_color_manual(values=c( "orange3","dodgerblue3"), name="Life history", labels=c("Annual", "Prennial"))+  scale_y_continuous(breaks = c(0,1), labels=c("       a","       p"))+
  geom_line(aes(x=droughtfreq, y=curves), col="black")+
  theme(axis.text.x = element_text(angle=90), strip.background = element_rect(fill="white"))+
labs(y="Life history", x="Drought frequency")+
  geom_text(data=fit_stats, aes(x=.4, y=.75, label=label), col="black", cex=5)

plot_grid(lines, fitplots, dates, ncol=1, labels="auto", rel_heights = c(1,1.2,1.5))

annuals<-subset(GBIF_cleaned, life_history=="a")
perennials<-subset(GBIF_cleaned, life_history=="p")

kout<-ks.test(annuals$days,perennials$days)
tout<-t.test(days~life_history, GBIF_cleaned)
bout<-bartlett.test(days~life_history, GBIF_cleaned)

dev.off()


## ----lineplots, fig.cap='(ref:lineplots)', out.width = "\\textwidth", fig.pos = "!h"----
knitr::include_graphics("../figures/line_and_dates.pdf", dpi = 108)

## ----create_r-references-------------------------------------------------
#r_refs(file = "r-references.bib")

## ----anovatable, results='asis', echo = F--------------------------------
apa_table(anovatable
          , caption = "Analysis of variance (ANOVA) to compare drought frequency as a function of life history, season, and their interaction while including species as a random effect.",
          digits=4
)

## ----speciesmeanstable, results = 'asis', echo = F-----------------------
species_means_table<-species_means
row.names(species_means_table) <- NULL
species_means_table<-species_means_table %>% dplyr::select(-life_history_num)
colnames(species_means_table)<-c("Species","LH","n","Winter","Spring","Summer","Fall")
apa_table(
  species_means_table
  , caption = "Heliophila species records and the mean drought frequencies during different seasons at the location of records "
  , note = "LH = Life history (a = annual, p = perennial). n=sample size of GBIF records. Seasons are mean drought frequencies observed at locations of records.",
  digits=2,
  font_size="small",
  longtable = TRUE
)

write.csv(species_means_table, "../data/species_means_table.csv")

## ----results = 'hide', echo = F------------------------------------------
write.csv(species_means_table, "../data/species_means_table.csv")

## ----make_drought_examples, results = 'hide', echo = F, warnings="hide"----
pdf("../figures/maps_drought_examples.pdf", width=3, height=6)

world <- map_data("world") # we already did this, but we can do it again
col_scale<-c("white","black")

drought<-raster::raster("~/VHI/raster/VHP.G16.C07.npp.P2015050.VH.nc.grd")
drought <-crop(drought, rbind(c(12,37),c(-35, -18)))

summer2015<-gplot(drought) + 
  geom_tile(aes(fill=value<4000)) + 
  scale_fill_manual(values = col_scale,  na.translate = F, name="VHI < 40") +
  geom_polygon(data = world, aes(x=long, y = lat, group=group),  fill = NA, color = "black") + 
  coord_fixed(xlim=c(13,35), ylim=c(-35, -22)) + 
  theme_classic()+
  theme(legend.position = c(0.92,.3),legend.text=element_text(size=rel(0.9)), legend.key.size=unit(.2,"line"), legend.background = element_rect(fill=alpha('gray90'), color = "black"), text = element_text(size=4.5))+
  labs(title="VHI detected drought during summer 2015", x="Longitude (°)", y="Latitude (°)")

nondrought<-raster::raster("~/VHI/raster/VHP.G16.C07.NJ.P1999050.VH.hdf.h5.grd")
nondrought <-crop(nondrought, rbind(c(12,37),c(-35, -18)))

summer1999<-gplot(nondrought) + 
  geom_tile(aes(fill=value<4000)) + 
  scale_fill_manual(values = col_scale,  na.translate = F, name="VHI < 40") +
  geom_polygon(data = world, aes(x=long, y = lat, group=group),  fill = NA, color = "black") + 
  coord_fixed(xlim=c(13,35), ylim=c(-35, -22)) + 
  theme_classic()+
  theme(legend.position = c(0.92,.3),legend.text=element_text(size=rel(0.9)), legend.key.size=unit(.2,"line"), legend.background = element_rect(fill=alpha('gray90'), color = "black"), text = element_text(size=4.5))+
  labs(title="VHI detected drought during summer 1999", x="Longitude (°)", y="Latitude (°)")

typical<-raster::raster("~/VHI/raster/VHP.G16.C07.NP.P2012050.VH.hdf.h5.grd")
typical <-crop(typical, rbind(c(12,37),c(-35, -18)))

summer2012<-gplot(typical) + 
  geom_tile(aes(fill=value<4000)) + 
  scale_fill_manual(values = col_scale,  na.translate = F, name="VHI < 40") +
  geom_polygon(data = world, aes(x=long, y = lat, group=group),  fill = NA, color = "black") + 
  coord_fixed(xlim=c(13,35), ylim=c(-35, -22)) + 
  theme_classic()+
  theme(legend.position = c(0.92,.3),legend.text=element_text(size=rel(0.9)), legend.key.size=unit(.2,"line"), legend.background = element_rect(fill=alpha('gray90'), color = "black"), text = element_text(size=4.5))+
  labs(title="VHI detected drought during summer 2012", x="Longitude (°)", y="Latitude (°)")

plot_grid(summer2015, 
          summer1999,
          summer2012, 
          labels = "auto", ncol = 1)


dev.off()


## ----mapsdroughtexamples, fig.cap='(ref:mapsdroughtexamples)'------------
knitr::include_graphics("../figures/maps_drought_examples.pdf", dpi = 108)

## ----make_spceiesmaps, results = 'hide', echo = F, warnings="hide"-------
pdf("../figures/speciesmaps.pdf", height=8, width=8)

world <- map_data("world") # we already did this, but we can do it again

plot_grid(plotlist=lapply(unique(GBIF_cleaned$species), function(x){
  
  lh<-unique((GBIF_cleaned %>% filter(species==x))$life_history)
  lh_col=ifelse(lh=="a", "orange3", "dodgerblue3")
  
  ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group=group),  fill = NA, color = "black") + 
    coord_fixed(xlim=c(13,35), ylim=c(-35, -22)) + 
    geom_point(GBIF_cleaned %>% filter(species==x), mapping=aes(x=decimalLongitude, y=decimalLatitude), pch=3, col=lh_col, alpha=0.3) +
    theme_classic()+
    theme( text = element_text(size=5))+
    labs(title=x, x=NULL, y=NULL)
  
}), ncol=6)

dev.off()


## ----speciesmaps, fig.cap='(ref:speciesmaps)'----------------------------
knitr::include_graphics("../figures/speciesmaps.pdf", dpi = 108)

