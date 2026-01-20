

# Function to run R script
run_r_script <- function(run_type, script_name, expression_file, ortho_dir, dups_file, exp_dir, add_pseudofunc, missing_expr_is_pseudo, exp_cutoff, PC, min_dups_per_species_pair, use_absolute_exp) {
  source(paste0("C:/Users/17735/Downloads/DuplicA/app/Scripts/", script_name))
  
  print(paste('is:', expression_file, ortho_dir, as.logical(add_pseudofunc), as.logical(missing_expr_is_pseudo), as.numeric(exp_cutoff), as.logical(PC), as.numeric(min_dups_per_species_pair), as.logical(use_absolute_exp)))
  print(paste('should be:', 'C:/Users/17735/Downloads/AAAAA___EXAMPLE_Expression.tsv', 'C:/Users/17735/Downloads/AAAAA_Results_Jan01', 'True', 'False', '1', 'False', '10', 'False'))
  
  tryCatch({
    if (run_type == 'OF') {
      main_CDROM(expression_file, ortho_dir, as.logical(add_pseudofunc), as.logical(missing_expr_is_pseudo), as.numeric(exp_cutoff), as.logical(PC), as.numeric(min_dups_per_species_pair), as.logical(use_absolute_exp))
    } else if (run_type == 'custom') {
      # does not work with custom input 
      #main_CDROM(dups_file, exp_dir, as.logical(add_pseudofunc), as.logical(missing_expr_is_pseudo), as.numeric(exp_cutoff), as.logical(PC), as.numeric(min_dups_per_species_pair), as.logical(use_absolute_exp))
    }
  }, error = function(e) {
    return(paste("An error occurred:\n", e$message))
  })
}


# allow input of directories 
dirInput <- function(inputId, label, width = NULL, buttonLabel = "Browse...", 
                     placeholder = "No directory selected", directory = TRUE) {
  
  restoredValue <- restoreInput(id = inputId, default = NULL)
  
  # Ensure the restored value is either NULL or a data frame.
  if (!is.null(restoredValue) && !is.data.frame(restoredValue)) {
    warning("Restored value for ", inputId, " has incorrect format.")
    restoredValue <- NULL
  }
  
  if (!is.null(restoredValue)) {
    restoredValue <- toJSON(restoredValue, strict_atomic = FALSE)
  }
  
  # Input element to allow directory selection (directory=TRUE is the key)
  inputTag <- tags$input(
    id = inputId,
    class = "shiny-input-file",
    name = inputId,
    type = "file",
    style = "position: absolute !important; top: -99999px !important; left: -99999px !important;",
    `data-restore` = restoredValue
  )
  
  # Allow directory selection if the `directory` argument is TRUE
  if (directory) {
    inputTag$attribs$webkitdirectory <- "webkitdirectory"
  }
  
  # Return the UI for file input that is actually selecting a directory
  div(class = "form-group shiny-input-container",
      style = css(width = validateCssUnit(width)),
      
      div(class = "input-group",
          tags$label(class = "input-group-btn input-group-prepend",
                     span(class = "btn btn-default btn-file",
                          buttonLabel,
                          inputTag
                     )
          ),
          tags$input(type = "text", class = "form-control",
                     placeholder = placeholder, readonly = "readonly"
          )
      ),
      
      tags$div(
        id=paste(inputId, "_progress", sep=""),
        class="progress active shiny-file-input-progress",
        tags$div(class="progress-bar")
      )
  )
}


# add CSS to style buttons make the sidebar narrower
add_css_style <- function() {
  
  tags$style(HTML("
    .sidebar-buttons .action-button {
      display: block;
      width: 100%;
      margin-bottom: 15px;
      background-color: #007BFF;  /* Blue buttons */
      border-color: #007BFF;
    }
    .sidebar-buttons .action-button:hover {
      background-color: #0056b3;
    }
  "))
}


# Page functions for each model - return content only, not full layouts
workflow_page <- function() {
  tagList(
    h3("Workflow Builder"),
    p("The workflow feature is available in the modern web interface."),
    p("To access the full interactive workflow builder with visual node-based construction:"),
    tags$ol(
      tags$li("Open your browser and navigate to ", tags$code("http://localhost:8000/Workflow")),
      tags$li("Or use the Gatsby frontend if running in development mode")
    ),
    tags$div(style = "margin-top: 20px;"),
    p("The Shiny interface provides individual model access through the sidebar buttons.")
  )
}

orthofinder_page <- function() {
  tagList(
    h3("OrthoFinder"),
    p("Phylogenetic orthology inference for comparative genomics."),
    p("Configure and run OrthoFinder analysis on your protein sequences."),
    tags$div(style = "margin-top: 20px;"),
    p("This page is under construction.")
  )
}

dnds_page <- function() {
  tagList(
    h3("dN/dS Analysis"),
    p("Calculate the ratio of nonsynonymous to synonymous substitutions."),
    p("This analysis helps detect selection pressure on duplicate genes."),
    tags$div(style = "margin-top: 20px;"),
    p("This page is under construction.")
  )
}

segregating_duplicates_page <- function() {
  tagList(
    h3("Segregating Duplicates"),
    p("Analyze segregating duplicate genes in populations."),
    p("Uses CNVSelectR to test for selection on copy number variants."),
    tags$div(style = "margin-top: 20px;"),
    p("This page is under construction.")
  )
}

expression_shift_page <- function() {
  tagList(
    h3("EVE Expression Shift"),
    p("Phylogenetic ANOVA for expression-based lineage divergence."),
    p("Detect shifts in gene expression patterns across evolutionary lineages."),
    tags$div(style = "margin-top: 20px;"),
    p("This page is under construction.")
  )
}

diversity_divergence_page <- function() {
  tagList(
    h3("EVE Diversity/Divergence"),
    p("Phylogenetic ANOVA expression-based selection test."),
    p("Test for selection using expression diversity and divergence."),
    tags$div(style = "margin-top: 20px;"),
    p("This page is under construction.")
  )
}

blat_page <- function() {
  tagList(
    h3("BLAT Analysis"),
    p("BLAST-Like Alignment Tool for sequence comparison."),
    p("Fast sequence alignment for identifying duplicate genes."),
    tags$div(style = "margin-top: 20px;"),
    p("This page is under construction.")
  )
}

blast_page <- function() {
  tagList(
    h3("BLAST Analysis"),
    p("Basic Local Alignment Search Tool."),
    p("Identify duplicate genes through sequence similarity search."),
    tags$div(style = "margin-top: 20px;"),
    p("This page is under construction.")
  )
}

cdrom_orthofinder_tab <- function() {
  tabPanel(
    "OrthoFinder Input",
    h4("CDROM with OrthoFinder Data"),
    p("Classification of Duplicate gene Retention Mechanisms using OrthoFinder output."),
    tags$div(style = "margin-top: 20px;"),
    p("This tab is under construction.")
  )
}

cdrom_custom_tab <- function() {
  tabPanel(
    "Custom Input",
    h4("CDROM with Custom Data"),
    p("Classification of Duplicate gene Retention Mechanisms using custom duplicate gene data."),
    tags$div(style = "margin-top: 20px;"),
    p("This tab is under construction.")
  )
}

