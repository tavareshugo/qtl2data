####
# Get Brachypodium Bd21 x Bd3-1 RIL example phenotype
###

library(readtext) # to read word document
library(readr)
library(stringr)
library(purrr)

#### download data ####

# Download word doc file
file_url <- "https://doi.org/10.1371/journal.pone.0038333.s002"
file_temp <- paste0(tempfile(), ".doc")

download.file(file_url, file_temp, method = "wget")

# Read text from word document
pheno_text <- readtext(file_temp)$text %>%
  read_lines()

# Keep only lines with actual data
# excluding lines segregating for resistance
pheno_text <- pheno_text[5:87]


#### Parse data ####

# Split text by "|"
pheno <- str_split(pheno_text, "\\|") %>%
  # Remove whitespaces
  map(str_trim) %>%
  # Extract id and phenotype (2 per line)
  map_dfr(function(i){
    data.frame(id = i[c(2, 5)],
               phenotype = i[c(3, 6)],
               stringsAsFactors = FALSE)
  })


# save RDS file
write_rds(pheno, "../Cui2012_phenotype.rds", compress = "gz")
