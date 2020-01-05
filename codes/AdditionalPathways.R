### This is a code for additional iPathwayGuide results
### I downloaded Pathway figures and gene tables from
### both KEGG pathways and GO terms
### The gene tables have all zero rows, so remove them
### and trim the file names: GO figures do not have
### term names, they only have GO IDs.

### pathway file directory
pathway_dir <- "C:/Research/CUMC/Hans_RNASeq/data/pathways_of_interest/"

### directory of gene tables and figures
result_dir <- "C:/Research/CUMC/Hans_RNASeq/results/pathway/pathway_images_of_interest/"

### load library
if(!require(xlsx, quietly = TRUE)) {
  install.packages("xlsx")
  require(xlsx, quietly = TRUE)
}

### get result sub-directory names
result_sub_dirs <- list.files(path = result_dir)

### for each sub-directory change the files in it
for(dir in result_sub_dirs) {
  
  ### for pathwyasTable: remove all-zero rows and change the column names in tables
  ###                    change the file names
  ### for goTable: change the file names
  if(grepl("pathwaysTable", dir)) {
    
    ### get csv files only
    files <- list.files(path = paste0(result_dir, dir), pattern = "\\.csv$")
    
    ### get core title only
    titles <- sapply(files, function(x) {
      temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
      temp <- paste(temp[-c(1:3)], collapse = "_")
      return(substr(temp, 1, nchar(temp)-4))
    })
    
    ### remove all-zero rows and change the column names in tables
    for(file in files) {
      
      ### load gene table
      gene_table <- read.table(file = paste0(result_dir, dir, "/", file),
                               header = TRUE, sep = ",",
                               stringsAsFactors = FALSE, check.names = FALSE)
      
      ### find all zero rows
      zero_rows <- sapply(1:nrow(gene_table), function(x) {
        return(gene_table[x,3] == 0 && gene_table[x,4] == 0 && gene_table[x,5] == 0)
      })
      
      ### remove all zero rows
      if(length(which(zero_rows)) > 0) {
        gene_table <- gene_table[!zero_rows,]
      }
      
      ### change the column names
      colnames(gene_table)[4:5] <- c("accumulation", "perturbation")
      
      ### save the gene table
      write.xlsx2(gene_table, file = paste0(result_dir, dir, "/", titles[file], ".xlsx"),
                  sheetName = titles[file], row.names = FALSE, col.names = TRUE)
      
      ### remove the csv file
      file.remove(paste0(result_dir, dir, "/", file))
      
    }
    
    ### get png files only
    files <- list.files(path = paste0(result_dir, dir), pattern = "\\.png$")
    
    ### change the names of the figure files
    for(file in files) {
      file.rename(from = paste0(result_dir, dir, "/", file),
                  to = paste0(result_dir, dir, "/", strsplit(file, split = "__", fixed = TRUE)[[1]][2]))
    }
    
  } else if(grepl("goTable", dir)) {
    
    ### load the GO result file
    result_file <- read.xlsx2(file = paste0(pathway_dir, dir, ".xlsx"),
                              sheetIndex = 1, row.names = 1)
    
    ### get file names
    files <- list.files(path = paste0(result_dir, dir))
    
    ### annotate the GO terms to the file names
    for(file in files) {
      
      ### get GO ID for the file
      go_id <- paste0("GO:", strsplit(file, split = "_", fixed = TRUE)[[1]][5])
      
      ### get GO term name for the ID
      go_term <- as.character(result_file[go_id, "GO_Name"])
      go_term <- paste(strsplit(go_term, split = " ", fixed = TRUE)[[1]], collapse = "_")
      
      ### change the file name
      if(grepl(".csv", file)) {
        ### load the file
        csv_file <- read.table(file = paste0(result_dir, dir, "/", file), header = TRUE,
                               stringsAsFactors = FALSE, check.names = FALSE, sep = ",")
        
        ### save the file with changed name
        write.xlsx2(csv_file, file = paste0(result_dir, dir, "/", go_term, "_geneTable.xlsx"),
                    sheetName = go_term, row.names = FALSE, col.names = TRUE)
        
        ### remove the csv file
        file.remove(paste0(result_dir, dir, "/", file))
      } else if(grepl(".png", file)) {
        ### change the name of the figure file
        file.rename(from = paste0(result_dir, dir, "/", file),
                    to = paste0(result_dir, dir, "/", go_term, "_goAncestry.png"))
      }
      
    }
    
  } else {
    stop("ERROR: the sub-directory should contain either pathwaysTable or goTable.")
  }
  
}
