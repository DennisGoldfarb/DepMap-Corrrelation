library(tidyverse)
library(here)

#-------------------------------------------------------------------------------
# Usage Rscript pairwiseCor.R [first_letter_of_gene]
#
# This calculates the pearson's correlation between each gene that starts with the
# specified first letter and all genes. It's split by first letter so you can
# parallelize this.
#-------------------------------------------------------------------------------


# Helper function. Use it to create a folder for each letter of the alphabet
createDir <- function(path)
{
  if (!file.exists(path))
  {
    dir.create(path)
  }
}

# Get the required first letter of the gene name
# This makes it easy to do it parallel on a compute cluster
# Which is nice because this takes a while
args <- commandArgs(trailingOnly=TRUE)
first_letter <- args[1]

# output path
base_path <- "Depmap CRISPR KO correlations"
letter_path <- str_c(base_path, "/", first_letter)
createDir(here(base_path))
createDir(here(letter_path))

# DepMap data
data_input <- read_csv(here("CRISPR_gene_effect.csv"), guess_max=2) %>%
  select(-DepMap_ID)

# Create empty output table
base_table <- tibble(gene=colnames(data_input)) %>%
  extract(gene, c("Gene Name", "Gene ID"), "(.*) \\((.*)\\)")



# For each gene, compute correlation with all other genes
for (i in seq(1, ncol(data_input)))
{
  gene1 <- colnames(data_input)[i]
  # filtered for genes that start with the required first letter
  if (str_sub(gene1,1,1) != first_letter) next

  print(gene1)

  gene_correlations <- base_table %>% mutate(pcc = 0)

  for (j in 1:ncol(data_input))
  {
    gene2 <- colnames(data_input)[j]
    pearsons <- cor(data_input[[gene1]], data_input[[gene2]])
    gene_correlations[j, 3] <- pearsons
  }

  write_csv(gene_correlations, here(str_c(letter_path,"/",gene1,".csv")))
}





