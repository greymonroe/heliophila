# Repository for developing the manuscript looking at life history in Heliophila in relation to drought frequency.

/data - contains raw input data files  
/figures contains figures created by scripts that will be included in the publication
/manuscript contains markdown, bibtex, and formating files to build manuscript
/pictures has a few pictures of specimens

to build the manuscript run
`knitr("manuscript/manuscript.rmd")`

you might have to install some packages first...
`
library("papaja");
library(raster);
library(tidyverse);
library(ggplot2);
library(logistf);
library(geiger) ;
library(ape) ;
library(phylolm);
library(gridExtra);
library(cowplot);
library(rasterVis);
library(ggtree);
library(grid);
library(png)
`



