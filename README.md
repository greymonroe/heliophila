# heliophila

/data - contains raw input data files  
/data/0006178-160311141623029/occurrence.csv - raw GBIF data

/src contains scripts to process and analyze raw data
/Rdata contains Rdata objects produced by scripts. 
/figures contains figures created by scripts that will be included in the publication
/tables contains tables created by scripts that will be included in the publication
/manuscript contains markdown, bibtex, and formating files to build manuscript

required packages
devtools::install_github("crsh/papaja")
install.packages("citr")

to build the manuscript all you need to do is run
knitr("manuscript/manuscript.rmd")



