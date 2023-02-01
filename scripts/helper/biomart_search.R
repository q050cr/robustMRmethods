

# get gene names from rs SNP ids

library(biomaRt)
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
# https://bio340files.s3.amazonaws.com/SearchBiomart.html  https://www.youtube.com/watch?v=bMkve2uBMZo&ab_channel=profbiot

# list databases (look at pulldown menu on website: choose database)
myMarts <- listMarts()

# biomart                version
# 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 108
# 2   ENSEMBL_MART_MOUSE      Mouse strains 108
# 3     ENSEMBL_MART_SNP  Ensembl Variation 108
# 4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 108

listEnsembl()   # the same since listMarts connects to Ensembl
# biomart                version
# 1         genes      Ensembl Genes 108
# 2 mouse_strains      Mouse strains 108
# 3          snps  Ensembl Variation 108
# 4    regulation Ensembl Regulation 108

# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# select database
ensembl <- useMart(biomart = "ENSEMBL_MART_SNP")
myDatasets <- listDatasets(ensembl)

# select dataset: hsapiens
ensembl <- useMart(biomart = "ENSEMBL_MART_SNP", 
                   dataset = "hsapiens_snp") #  defaults: host = "https://www.ensembl.org"; path="/biomart/martservice"
                   
# List filters
myFilters <- listFilters(ensembl)

# make a vector for your filters (here all rs_ids)
filter1 <- c('snp_filter')

# make a vector for the values that will go into your filter
values1 <- my_data$SNP

# list attributes (columns we download from biomart)
myAttributes <- listAttributes(ensembl)
# make a variable for your attributes
att1 <- c('refsnp_id', 'minor_allele', 'minor_allele_freq', 
          'phenotype_name', 
          'phenotype_description'  # several descriptions per SNP possible!
          )

att2 <- c('refsnp_id', 'minor_allele', 'minor_allele_freq', 
          'clinical_significance', 'associated_gene', 
          'phenotype_name', 
          'ensembl_gene_stable_id', 'ensembl_gene_name'
          )

# SEARCH ------------------------------------------------------------------
# search results

searchResults1 <- getBM(attributes = att1, 
                       filters = filter1,
                       values = values1, mart = ensembl)

searchResults2 <- getBM(attributes = att2, 
                        filters = filter1,
                        values = values1, mart = ensembl)

saveRDS(searchResults2, "data/biomartBMIsearch.rds")
