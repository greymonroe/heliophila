
pdf("../figures/phylogeny.pdf", height=6, width=5)

treeplot<-ggtree(tree, layout = "rectangular")+
  xlim(c(0,2))+
  geom_tippoint(size=3,col=ifelse(species_life_history$life_history[match(tree$tip.label,species_life_history$species)]=="a", "orange3", "dodgerblue3"))+
  geom_tiplab(offset = .01 )

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
                labels = c("b","c","d","e"), ncol = 1)

plot_grid(treeplot, pics, ncol=2,labels=c("a",""), rel_widths = c(1,0.35))


dev.off()
