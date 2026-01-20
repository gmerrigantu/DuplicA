
get_data_for_organisms <- function(selected_organisms, data_types, selected_database_protein, selected_database_cds, selected_database_genome, prot_data_dir, cds_data_dir, genome_data_dir, must_be_reference) {
  
  if ('Proteomes' %in% data_types){
    getProteomeSet(db = selected_database_protein, organism = selected_organisms, path = prot_data_dir, reference = must_be_reference)
  }

  if ('CDS' %in% data_types){
    getCDSSet(db = selected_database_cds, organism = selected_organisms, path = cds_data_dir, reference = must_be_reference)
  }
  
  if ('Genomes' %in% data_types){
    getGenomeSet(db = selected_database_genome, organism = selected_organisms, path = genome_data_dir, reference = must_be_reference)
  }
  
  
}


get_organisms_prot_fasta_data <- function(selected_organism, prot_files, selected_database_protein) {
  if (selected_database_protein == 'refseq') {
    prot_file <- prot_files[
      grepl(selected_organism, prot_files) &
        grepl("_protein_refseq\\.faa(\\.gz)?$", prot_files)
    ]
  } else {
    prot_file <- prot_files[grepl(selected_organism, prot_files) & grepl("\\.pep\\.", prot_files)]
  }
  prot_fasta_data <- readAAStringSet(prot_file)
  
  return(prot_fasta_data)
}


get_organisms_cds_fasta_data <- function(selected_organism, cds_files, selected_database_cds) {
  if (selected_database_cds == 'refseq') {
    cds_file <- cds_files[
      grepl(selected_organism, cds_files) &
        grepl("_cds_from_genomic\\.fna(\\.gz)?$", cds_files)
    ]
  } else {
    cds_file <- cds_files[grepl(selected_organism, cds_files) & grepl("\\.cds\\.", cds_files)]
  }
  cds_fasta_data <- readDNAStringSet(cds_file)
  
  return(cds_fasta_data)
  
}


get_prot_transcript_seq <- function(prot_fasta_data, keep_which_transcript, selected_organism) {
  prot_seqs_df <- as.data.frame(prot_fasta_data) %>%
    rownames_to_column('names') %>%
    mutate(
      protein_id = str_extract(names, "^[^ ]+"),
      gene_id = str_extract(names, "(?<=gene:)[^ ]+"),
      transcript_id = str_extract(names, "(?<=transcript:)[^ ]+")
    ) %>%
    dplyr::select(gene_id, transcript_id, protein_id, prot_seq = x) %>%
    mutate(len = nchar(prot_seq)) 
  
  
  any(duplicated(prot_seqs_df$transcript_id)) # should be FALSE
  any(duplicated(prot_seqs_df$protein_id)) # should be FALSE
  
  if (keep_which_transcript == 'longest'){
    prot_seqs_df <- prot_seqs_df %>%
      group_by(gene_id) %>%
      filter(len == max(len)) %>% 
      slice_head(n = 1) %>%  # keep only one even if max length is the same 
      ungroup() 
  }
  
  if (keep_which_transcript == 'first'){
    prot_seqs_df <- prot_seqs_df %>%
      group_by(gene_id) %>%
      slice_head(n = 1) %>%  
      ungroup() 
  }
  any(duplicated(prot_seqs_df$gene_id)) # should be FALSE
  
  gene_transcript <- prot_seqs_df %>% dplyr::select(gene_id, transcript_id)
  selected_organism <- gsub(' ', '_', selected_organism)
  write.table(gene_transcript, file = paste0(here_results, '/Fastas/kept_transcript/', selected_organism, '_transcript_kept_per_gene.tsv'), sep = '\t', quote = F, row.names = F)
  
  return(prot_seqs_df)
}


add_cds_transcript_seq <- function(cds_fasta_data, prot_seqs_df) {
  
  cds_seqs_df <- as.data.frame(cds_fasta_data) %>%
    rownames_to_column('names') %>%
    mutate(transcript_id = str_extract(names, "^[^ ]+")) %>%
    dplyr::select(transcript_id, cds_seq = x)
  
  seqs_df <- merge(cds_seqs_df, prot_seqs_df, by = 'transcript_id')
  
  seqs_df <- seqs_df %>% dplyr::select(gene_id, transcript_id, protein_id, cds_seq, prot_seq)
  
  return(seqs_df)
}


create_output_dirs <- function() {
  prot_output_dir <- paste0(here_results, '/Fastas/Protein_Fastas/')
  dir.create(prot_output_dir, recursive = T)
  cds_output_dir <- paste0(here_results, '/Fastas/Nucleotide_Fastas/')
  dir.create(cds_output_dir, recursive = T)
  
  return(list(prot_output_dir = prot_output_dir, cds_output_dir = cds_output_dir))
}


seq_df_to_fasta_files <- function(seqs_df, prot_output_dir, cds_output_dir, selected_organism) {
  
  prot_output_file <- paste0(prot_output_dir, selected_organism, '_prot.fasta')
  cds_output_file <- paste0(cds_output_dir, selected_organism, '_cds.fasta')
  
  seqs_df %>%
    mutate(fasta = paste0('>', gene_id, '\n', prot_seq)) %>%
    pull(fasta) %>%
    writeLines(con = prot_output_file)
  
  seqs_df %>%
    mutate(fasta = paste0('>', gene_id, '\n', cds_seq)) %>%
    pull(fasta) %>%
    writeLines(con = cds_output_file)
}


move_genome_files <- function(selected_organism, genome_files, selected_database_genome) {
  
  
  # get the genome of the given species 
  if (selected_database_genome == 'refseq') {
    genome_file <- genome_files[
      grepl(selected_organism, genome_files) &
        grepl("_genomic\\.fna(\\.gz)?$", genome_files)
    ]
  } else {
    genome_file <- genome_files[grepl(selected_organism, genome_files) & grepl("\\.dna\\.", genome_files)]
  }
  
  if (length(genome_file) == 0) {
    print(paste0('No genome file found for ', selected_organism))
    return()
  }
  
  # create a new directory for the genome files
  genome_output_dir <- paste0(here_results, '/Fastas/Genome_Fastas/')
  dir.create(genome_output_dir, recursive = T)
  
  # copy the file to the output directory
  file.copy(genome_file, new_genome_file)
  
  # rename the new file 
  new_genome_file <- paste0(genome_output_dir, selected_organism, '_genome.fasta')
  
  
  print(genome_file)
  print(length(genome_file))
  print(new_genome_file)
  print(length(new_genome_file))

  file.rename(genome_file, new_genome_file)
  file.remove(paste0(genome_output_dir, '/', basename(genome_file)))
  
  
  
}


main_public_datasets <- function(selected_organisms, data_types, selected_database_protein, selected_database_cds, selected_database_genome, keep_which_transcript, must_be_reference) {
  
  if (missing(selected_database_protein) || is.null(selected_database_protein)) {
    selected_database_protein <- 'refseq'
  }
  if (missing(selected_database_cds) || is.null(selected_database_cds)) {
    selected_database_cds <- 'refseq'
  }
  if (missing(selected_database_genome) || is.null(selected_database_genome)) {
    selected_database_genome <- 'refseq'
  }
  
  prot_data_dir <- paste0(here_results, '/public_datasets_output/protein_data')
  cds_data_dir <- paste0(here_results, '/public_datasets_output/cds_data')
  genome_data_dir <- paste0(here_results, '/public_datasets_output/genome_data')
  
  # get the data types from the public datasets 
  get_data_for_organisms(selected_organisms, data_types, selected_database_protein, selected_database_cds, selected_database_genome, prot_data_dir, cds_data_dir, genome_data_dir, must_be_reference)
  
  # list all files provided by the public datasets 
  prot_files <- list.files(prot_data_dir, full.names = TRUE)
  cds_files <- list.files(cds_data_dir, full.names = TRUE)
  genome_files <- list.files(genome_data_dir, full.names = TRUE)

  selected_organisms <- gsub(' ', '_', selected_organisms)
  
  # interate over each selected organism 
  for (selected_organism in selected_organisms) {
    
    # read in the fasta files for the given organism 
    prot_fasta_data <- get_organisms_prot_fasta_data(selected_organism, prot_files, selected_database_protein)
    cds_fasta_data <- get_organisms_cds_fasta_data(selected_organism, cds_files, selected_database_cds)
    
    # format the fasta data into a table. keep only one sequence per gene
    prot_seqs_df <- get_prot_transcript_seq(prot_fasta_data, keep_which_transcript = 'longest', selected_organism)
    seqs_df <- add_cds_transcript_seq(cds_fasta_data, prot_seqs_df)
    
    # write the sequences into fasta files 
    dirs <- create_output_dirs()
    seq_df_to_fasta_files(seqs_df, dirs$prot_output_dir, dirs$cds_output_dir, selected_organism)
    
    
    # write genome fastas if selected
    if ('Genomes' %in% data_types) {
      move_genome_files(selected_organism, genome_files, selected_database_genome)
    }
    
  }
  
  # delete the public_datasets_output directory
  unlink(paste0(here_results, '/public_datasets_output'), recursive = T)
  
}





# deal with errors from when data can not be found
# deal with warnings caused by creating already existing folders
# allow all databases to be checked for data
