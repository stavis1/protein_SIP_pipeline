install.packages("BiocManager")
library("BiocManager")
BiocManager::install(c("biocgeneric", 
    "biostrings", 
    'iranges', 
    's4vectors', 
    'xvector', 
    'zlibbioc'),
    ask = FALSE)
