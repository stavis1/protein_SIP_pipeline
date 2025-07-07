if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.8")
BiocManager::install(c("biocgeneric", 
    "biostrings", 
    'iranges', 
    's4vectors', 
    'xvector', 
    'zlibbioc'))

install.packages(c("cli", 
    'crayon', 
    'dplyr', 
    'ellipsis', 
    'fansi', 
    'generics', 
    'glue', 
    'lifecycle', 
    'magrittr', 
    'pillar', 
    'pkgconfig', 
    'purrr', 
    'r6', 
    'rlang', 
    'stringi', 
    'stringr',
    'tibble',
    'tidyr',
    'tidyselect',
    'utf8',
    'vctrs'))
