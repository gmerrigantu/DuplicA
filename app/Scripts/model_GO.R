

#source('./app/Scripts/multispecies_functions.R')


library(biomaRt)
# Load clusterProfiler - optional to avoid startup errors
tryCatch({
  library(clusterProfiler)
}, error = function(e) {
  # Package not available - will fail when GO functions are called
  NULL
})
library(biomartr)




main_go <- function(all_copies, file_organism_table) {
  colnames(file_organism_table)[2] <- 'organism_scientific_name'
  file_organism_table$organism_scientific_name <- gsub('_', ' ', file_organism_table$organism_scientific_name)
  
  
  # iterate over all organisms with data 
  all_go_output <- data.frame()
  for (chosen_organism in file_organism_table$organism_scientific_name) {
    
    # get file name for the chosen organism
    chosen_protein_file_name <- get_protein_file_name(chosen_organism, file_organism_table)
    
    # get all genes for the organism (any duplicate copies and any ancestral copies)
    genes <- get_genes_for_organism(chosen_protein_file_name)
    
    # check if data is available for the given organism
    avail_data <- get_avail_data_for_organism(chosen_organism, topic = 'go_id')
    if (nrow(avail_data) == 0) {return(paste0('No GO data is available for '), chosen_organism)}
    
    
    # get gene ontology data for the organism, when only one dataset is available
    if (nrow(avail_data) == 1) {
      go_output <- getGO(organism = chosen_organism, 
                         genes = genes$gene, 
                         filters = 'ensembl_gene_id')
      
      go_output <- go_output %>% 
        select(gene_id = ensembl_gene_id, 
               go_id = goslim_goa_accession,
               go_description = goslim_goa_description)
      
    }
    
    # use the first dataset when multiple are available 
    if (nrow(avail_data) > 1) {
      chosen_mart <- avail_data$mart[1] # CHANGEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEee
      chosen_dataset <- avail_data$dataset[1] 
      
      chosen_data <- useDataset(dataset = chosen_dataset, 
                                mart = useMart(chosen_mart))
      
      go_output <- getBM(attributes = c('ensembl_gene_id', "go_id", 'namespace_1003', 'goslim_goa_description', 'name_1006'),
                         filters = "ensembl_gene_id",
                         values = genes$gene,
                         mart = chosen_data)
      
      go_output <- go_output %>% 
        select(gene_id = ensembl_gene_id, 
               go_id = go_id,
               go_description = goslim_goa_description) 
      
      
      
      
    }
    
    # combine go output for all species 
    go_output <- go_output %>% mutate(protein_file_name = chosen_organism)
    all_go_output <- rbind(all_go_output, go_output)
  }
  
  all_copy_go <- left_join(all_copies, all_go_output, by = c('gene' = 'gene_id', 'protein_file_name'))
  
  go_output_path <- paste0(here_results, '/Gene_Ontology.tsv')
  write.table(all_copy_go, file = go_output_path, row.names = F, sep = '/t')
}





compare_copy_go <- function(all_copy_go) {
  
  
  t <- all_copy_go %>%
    group_by(gene, go_id) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Orthogroup) %>%
    filter(!any(is.na(across(everything())))) %>%
    summarize(
      all_d1_d2 = all(go_id[copy == "dup_1"] %in% go_id[copy == "dup_2"]),
      all_d2_d1 = all(go_id[copy == "dup_2"] %in% go_id[copy == "dup_1"]),
      all_d1_anc = all(go_id[copy == "dup_1"] %in% go_id[copy == "ancestral_copy"]),
      all_anc_d1 = all(go_id[copy == "ancestral_copy"] %in% go_id[copy == "dup_1"]),
      all_d2_anc = all(go_id[copy == "dup_2"] %in% go_id[copy == "ancestral_copy"]),
      all_anc_d2 = all(go_id[copy == "ancestral_copy"] %in% go_id[copy == "dup_2"]),
      any_d1_d2 = any(go_id[copy == "dup_1"] %in% go_id[copy == "dup_2"]),
      any_d2_d1 = any(go_id[copy == "dup_2"] %in% go_id[copy == "dup_1"]),
      any_d1_anc = any(go_id[copy == "dup_1"] %in% go_id[copy == "ancestral_copy"])
    )
  
}





# find genes associated with specific go 

# System-biology level: GO.db

# Genome centric GenomicFeatures packages include
# Transcriptome level: e.g. TxDb.Hsapiens.UCSC.hg19.knownGene, EnsDb.Hsapiens.v75

