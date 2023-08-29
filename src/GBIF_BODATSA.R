library(raster)
library(tidyverse)
library(data.table)
library(ggplot2)

GBIF_BODATSA<-fread("data/GBIF_GRIDREF.csv")


cor.test(GBIF_BODATSA$lat.x, GBIF_BODATSA$lat.y)
cor.test(GBIF_BODATSA$lon.x, GBIF_BODATSA$lon.y)

pdf("figures/lat_lon_comparison.pdf", width=1.5, height=1.5)
ggplot(all_means, aes(x=lat.x, y=lat.y))+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="Latitude (GBIF)")+
  scale_y_continuous(name="Latitude (BODATSA)")

ggplot(all_means, aes(x=lon.x, y=lon.y))+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="Longitude (GBIF)")+
  scale_y_continuous(name="Longitude (BODATSA)")
dev.off()

GBIF_BODATSA_long<-fread("data/GBIF_GRIDREF_long.csv")


pdf("figures/map_comparison.pdf", width=3, height=2.5)
world <- map_data("world") 
ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group=group),  fill = NA, color = "black") + 
  coord_fixed(xlim=c(13,35), ylim=c(-35, -22)) + 
  geom_point(all_means2, mapping=aes(x=lon, y=lat, fill=src, group=species), alpha=1, pch=21) +
  geom_line(all_means2, mapping=aes(x=lon, y=lat, group=species), alpha=1) +
  theme_classic()+
  theme( text = element_text(size=8), legend.position = "top")+
  labs(x="Longitude (°)", y="Latitude (°)")+
  scale_fill_manual(values = c("black","gray"), name="Source")
dev.off()


GBIF<-fread("data/GBIF.csv")
BODATSA<-fread("data/GRIDREF.csv")

all<-rbind(BODATSA, GBIF2)

all$loc<-paste(all$decimalLatitude, all$decimalLongitude, all$species)

all$obs<-paste(all$decimalLatitude, all$decimalLongitude, all$species)

BODATSAinGBIF<-sapply(unique(all[src=="GRIDREF"]$obs), function(u){
  u %in% all[src=="GBIF"]$obs
  
}) 
prop.table(table(BODATSAinGBIF))

GBIFinBODATSA<-sapply(unique(all[src=="GBIF"]$obs), function(u){
  u %in% all[src=="GRIDREF"]$obs
  
}) 
prop.table(table(GBIFinBODATSA))


# plot species means ------------------------------------------------------

life_history_drought_means<- species_means %>% 
  dplyr::select(life_history, Winter_mean, Spring_mean, Summer_mean,Autumn_mean)  %>%
  gather(variable, value, -life_history) %>%
  group_by(life_history, variable) %>%
  summarise(mean = mean(value), se = sd(value)/sqrt(length(value))) %>%
  mutate(variable = reorder(variable, c(4,2,3,1)))

lines<-ggplot(life_history_drought_means, aes(x=life_history, col=life_history,group=life_history, y=mean))+
  theme_bw(base_size = 6)+
  geom_point(size=2)+
  #geom_line(lwd=0.5, lty=1)+
  facet_grid(~variable, labeller = as_labeller(c("Winter_mean"="Winter","Spring_mean"="Spring","Summer_mean"="Summer","Autumn_mean"="Autumn")))+
  geom_errorbar(width=0,lwd=1, aes(ymin=mean-se, ymax=mean+se))+
  scale_y_continuous(limits=c(0.25,0.45))+
  scale_x_discrete(labels=NULL)+
  scale_color_manual(values=c( "orange3","dodgerblue3"), name="Life history", labels=c("Annual", "Prennial"))+
  labs(x=NULL,y="Drought Frequency")+
  theme(axis.ticks.x =element_blank(), strip.background =element_rect(fill="white"))

pdf("figures/BODATSA_means.pdf", width=4, height=2)
lines
dev.off()
# Phylogenetic comparison -------------------------------------------------

species_means <- BODATSA %>%
  dplyr::select(species, life_history, Winter, Summer, Spring, Autumn) %>%
  group_by(species, life_history) %>%
  summarise(n=n(), Winter_mean=mean(Winter), Spring_mean=mean(Spring), Summer_mean=mean(Summer), Autumn_mean=mean(Autumn)) %>% 
  as.data.frame()


tree<-read.nexus("data/TrimmedMCCwoPolytomies.nex")
tree<- rescale(tree, "depth", 1)
tree<-keep.tip(tree, tip = gsub(" ","_", as.character(species_means$species)))

species_means$life_history_num<-as.numeric(as.factor(species_means$life_history))-1
row.names(species_means)<-gsub(" ", "_", species_means$species)
seasons<-c("Winter_mean", "Spring_mean", "Summer_mean","Autumn_mean")

phylomodelstable<-data.frame()

for(s in 1:length(seasons)){
  fit1 <-  phyloglm(species_means$life_history_num~species_means[,seasons[s]], data=species_means, phy=tree)
  sum<-summary(fit1)
  row.names(sum$coefficients)<-NULL
  sum_df<-data.frame(Predictor=c("Intercept", gsub("_mean", " drought freq.", seasons[s])), sum$coefficients)
  phylomodelstable<-rbind(phylomodelstable, sum_df)
}

phylomodelstable <- phylomodelstable[,c( "Predictor","Estimate","p.value")]
phylomodelstable


# plot phylogenetic test --------------------------------------------------

fit <-  phyloglm(species_means$life_history_num~species_means$Spring_mean, data=species_means, phy=tree)
cc <- coef(fit)
spcurve<-plogis(cc[1]+cc[2]*sort(species_means$Spring_mean))

fit <-  phyloglm(species_means$life_history_num~species_means$Summer_mean, data=species_means, phy=tree)
cc <- coef(fit)
scurve<-plogis(cc[1]+cc[2]*sort(species_means$Summer_mean))

fit <-  phyloglm(species_means$life_history_num~species_means$Autumn_mean, data=species_means, phy=tree)
cc <- coef(fit)
fcurve<-plogis(cc[1]+cc[2]*sort(species_means$Autumn_mean))

fit <-  phyloglm(species_means$life_history_num~species_means$Winter_mean, data=species_means, phy=tree)
cc <- coef(fit)
wicurve<-plogis(cc[1]+cc[2]*sort(species_means$Winter_mean))

curvedataframe<-data.frame(
  droughtfreq=c(sort(species_means$Winter_mean),
                sort(species_means$Spring_mean),
                sort(species_means$Summer_mean),
                sort(species_means$Autumn_mean)),
  season=rep(c("Winter","Spring","Summer","Autumn"), each=42),
  curves=c(wicurve, spcurve, scurve, fcurve),
  lifehistory=c(species_means$life_history_num[order(species_means$Winter_mean)],
                species_means$life_history_num[order(species_means$Spring_mean)],
                species_means$life_history_num[order(species_means$Summer_mean)],
                species_means$life_history_num[order(species_means$Autumn_mean)])
)
curvedataframe$lifehistory<-curvedataframe$lifehistory
curvedataframe$season<-factor(curvedataframe$season, levels=c("Winter","Spring","Summer","Autumn"))

fit_stats<-data.frame(season=factor(c("Winter","Spring","Summer","Autumn"),levels=c("Winter","Spring","Summer","Autumn")), label=c("","*","**","**"), lifehistory="1")

fitplots<-ggplot(curvedataframe, aes(x=droughtfreq, col=as.factor(lifehistory), y=lifehistory))+
  geom_jitter(height = .03, width=0, alpha=0.7)+
  facet_grid(~season)+
  theme_bw(base_size = 6)+
  scale_color_manual(values=c( "orange3","dodgerblue3"), name="Life history", labels=c("Annual", "Prennial"))+  scale_y_continuous(breaks = c(0,1), labels=c("       A","       P"))+
  geom_line(aes(x=droughtfreq, y=curves), col="black")+
  theme(axis.text.x = element_text(angle=90), strip.background = element_rect(fill="white"))+
  labs(y="Life history", x="Drought frequency")+
  geom_text(data=fit_stats, aes(x=.4, y=.75, label=label), col="black", cex=5)

pdf("figures/BODATSA_fit.pdf", width=4, height=2)
fitplots
dev.off()
# ANOVA comparison --------------------------------------------------------

BODATSA$life_history<-as.factor(BODATSA$life_history)
melted<-BODATSA %>% select(life_history, Winter, Spring, Summer,Autumn, species) %>% gather(Winter, Spring, Summer,Autumn, key="season",value="droughtfreq")

fit<-lmer(droughtfreq~life_history*season+(1|species), melted)
anovatable<-data.frame(predictor=c("life history","season","life history x season"), anova(fit))
colnames(anovatable)[7]<-"p-value"
row.names(anovatable)<-NULL

post_hoc<-emmeans(fit, list(pairwise ~ life_history*season), adjust = "tukey")
post_hoc<-as.data.frame(post_hoc$`pairwise differences of life_history, season`)
post_hoc<-post_hoc[c(1,14,23,28), ]
post_hoc$contrast<-c("Autumn","Spring","Summer","Winter")

post_hoc$p.value

# boxplots ----------------------------------------------------------------

BODATSA_long<- BODATSA %>%
  dplyr::select(species, life_history, Winter, Spring, Summer,Autumn) %>%
  gather(season, drought_freq, Winter, Spring, Summer,Autumn)

BODATSA_long$season<-factor(BODATSA_long$season, levels=c("Winter","Spring","Summer","Autumn"))

boxplots<-ggplot(BODATSA_long, aes(x=life_history, col=life_history, y=drought_freq))+
  facet_wrap(~season, ncol=4 )+
  geom_jitter(width = .3, alpha=0.05, pch=3)+
  geom_boxplot(outlier.colour = "white", outlier.fill = NA, fill=NA, width=.3)+
  scale_color_manual(values=c("orange3","dodgerblue3"))+
  theme_classic(base_size = 6)+
  theme(legend.position = "none", text = element_text(size=8))+
  scale_x_discrete(name="Life history", labels=c("Annual","Perennial"))+
  labs(y="Drought frequency")

pdf("figures/BODATSA_boxplots.pdf", width=4, height=2)
boxplots
dev.off()



