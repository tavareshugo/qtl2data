####
# Get Brachypodium Bd21 x Bd3-1 RIL genotypes
###

library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(stringr)


#### Read data ####

# Download xls file
snp_file_url <- "https://doi.org/10.1371/journal.pone.0038333.s004"
snp_file_temp <- tempfile()

download.file(snp_file_url, snp_file_temp, method = "wget")

# Read genotype data (each chromosome is on a different sheet)
snp <- lapply(paste0("Chr", 1:5), function(chr){
  # In Chr1 sheet need to skip first 4 lines, otherwise 2 lines only
  if(chr == "Chr1"){
    skip_n <- 4
  } else {
    skip_n <- 2
  }

  read_excel(snp_file_temp,
             sheet = chr,
             skip = skip_n,
             col_types = "text")
})

# Bind list into single data.table
snp <- bind_rows(snp, .id = "chrom")

# Adjust column names
snp <- snp %>%
  rename(marker = Locus,
         pos_cm = `Genetic Position (cM)`,
         pos_bp = `Physical location of SNP position`,
         comments = Comments)

# Adjust column types
snp <- snp %>%
  mutate(pos_bp = as.numeric(pos_bp),
         pos_cm = as.numeric(pos_cm))

# Convert to long format
snp <- snp %>%
  gather("id", "genotype", matches("^RIL"))

# Assumming "-" is missing genotype
snp <- snp %>%
  mutate(genotype = ifelse(genotype == "-", NA, genotype))

# add columns regarding types of comments
snp <- snp %>%
  mutate(order_disagrees = ifelse(is.na(comments), FALSE, str_detect(comments, "order disagrees")),
         in_f2_map = ifelse(is.na(comments), TRUE, !str_detect(comments, "F2 map")),
         small_lg = ifelse(is.na(comments), FALSE, str_detect(comments, "small linkage group")))

# save compressed table of snps
write_rds(snp, "../Cui2012_genotypes.rds", compress = "gz")


#### Clean SNP dataset ####

# Remove Bsr1 marker (no physical position)
snp <- snp %>% filter(marker != "Bsr1")

# Remove BD3399_2 marker (typing same SNP as BD3399_1)
snp <- snp %>% filter(marker != "BD3399_2")

# Remove BD0419_4 marker (very high heterozygosity)
snp <- snp %>% filter(marker != "BD0419_4")


#### Create R/qtl2 objects ####

# Write genotype file
snp %>%
  select(id, marker, genotype) %>%
  spread(marker, genotype) %>%
  write_csv("../Cui2012_qtl2geno.csv")

# Write genetic map file
snp %>%
  distinct(marker, chrom, pos_cm) %>%
  write_csv("../Cui2012_qtl2gmap.csv")

# Write physical map file
snp %>%
  distinct(marker, chrom, pos_bp) %>%
  write_csv("../Cui2012_qtl2pmap.csv")

# Create control file - consider heterozygous sites as missing
file.remove("../Cui2012_qtl2cross.yaml")
qtl2::write_control_file(output_file = "../Cui2012_qtl2cross.yaml",
                         crosstype = "riself",
                         geno_file = "Cui2012_qtl2geno.csv",
                         gmap_file = "Cui2012_qtl2gmap.csv",
                         pmap_file = "Cui2012_qtl2pmap.csv",
                         geno_codes = c(a = 1, b = 2),
                         alleles = c("a", "b"),
                         na.strings = c("NA", "h"),
                         sep = ",",
                         description = "Brachypodium distachyon Bd3-1 x Bd21 RIL population, data from Cui et al. 2012")

# Zip data
setwd("..")
qtl2::zip_datafiles("Cui2012_qtl2cross.yaml", quiet = FALSE, overwrite = TRUE)

# remove uncompressed files
unlink("Cui2012_*.csv")
unlink("Cui2012_*.yaml")
