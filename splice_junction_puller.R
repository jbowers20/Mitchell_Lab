#################### User Input ####################
## dir where splice_junction_puller and Utils files are located ##

dir <- '~/Documents/splice_junction_project'

## Set output location
output_loc <- '~/Documents/splice_junction_project/analytic_docs'

## Set directory where input files are located ##
input_loc <- '~/Documents/splice_junction_project/input_files'
####################################################

## Finds all GFF3 files in specified directory ##
setwd(input_loc)
file_list <- list.files(pattern = ".gff")
## Imports ## 
library('spgs') ## Used to reverse compliment... you may need to install

## Loads Utils package ##
source(paste0(dir, '/Utils.R')) ## file must be in working directory

# For each GFF3 file located in this directory, do the following...##
for (file in file_list) {
  setwd(input_loc) 
  
  ## Reads the entire file ##
  file_lines <- readLines(file)
  fasta_line <- 0 ## Marks begining of fasta sequence
  fasta_seq <- "" ## Keeps track of fasta sequence
  comment_num <- 0 ## Keeps track of the number of comments in the file prior to the fasta section
  
  for (i in 1:length(file_lines)) {
    if (substr(file_lines[i], 1, 1) == ">") {
      fasta_line <- i
    }
    
    if (substr(file_lines[i], 1, 1) == "#") {
      comment_num <- comment_num + 1 
    }
    
    ## Read fasta seq ##
    if (fasta_line != 0 && i > fasta_line) {
      fasta_seq <- paste0(fasta_seq, file_lines[i])
    }
  }
  
  ## Create a data frame of annotations and pulls name from attributes section ##
  df <- gffRead(file, nrows = fasta_line - comment_num - 1)
  df <- df[df$feature == "exon", ] ## Removes all non-exon elements
  df <- df[!grepl(pattern = "Fragment", x= df$attributes, ignore.case = TRUE), ] ## Removes all exon fragments
  df <- df[!grepl(pattern = "Truncated", x= df$attributes, ignore.case = TRUE), ] ## Removes all truncated exons
  cat(' ...', dim(df)[1], "exons\n") ## prints number of exons
  df$Name <- gsub("%2C", ",", getAttributeField(df$attributes, "Name")) ## creates new row in df called Name and subs "%2C" for ","
  
  ## Create table to print
  output_df <- data.frame(gene = character(), exon = character(), 
                        acceptor = character(), donor = character())
  
  
  for (i in 1:dim(df)[1]) {
    
    if (dim(df)[1] != 0) { ## checks so that the size of the df is not equal to 0
      ## Define row ##
      row <- df[i,]
      
      ## Define gene and exon letter ##
      vec <- strsplit(row$Name, ", ")
      gene <- vec[[1]][1]
      exon <- vec[[1]][2]
      
      
      ## Find Donor (5') Site ##
      donor_site <- NA
      if (!grepl("E$", row$Name)) { ## No Donor Site for E
        if (row$strand == "+") {
          donor_site <- substr(fasta_seq, row$end + 1, row$end + 40)
        } else if (row$strand == "-") {
          donor_site <- reverseComplement(substr(fasta_seq, row$start - 40, row$start - 1), "dna", "as is")
        } else {
          stop("Error: Invalid strand character")
        }
      }
      
      
      ## Find Acceptor (3') Site ##
      acpt_site <- NA
      if (!grepl("A1$|A$|A1\\?$|A\\?$", row$Name)) { ## No Acceptor Site for A1/A or A1?/A?
        if (row$strand == "+") {
          acpt_site <- substr(fasta_seq, row$start - 40, row$start - 1)
        } else if (row$strand == "-") {
          acpt_site <- reverseComplement(substr(fasta_seq, row$end + 1, row$end + 40), "dna", "as is")
        } else {
          stop("Error: Invalid strand character")
        }
      }
      
      output_df[dim(output_df)[1] + 1,] <- c(gene,exon, acpt_site,donor_site) 
    }
    
  }
  ## Writes output to text file ##
  setwd(output_loc)
  if (dim(output_df)[1] > 0) {
    ## write.table(output_df, paste0(df[1,]$seqname,"_splice_sites.txt"), row.names = F, sep = "\t")
    saveRDS(output_df, paste0(df[1,]$seqname,"_splice_sites.rds"))
  } else {
    cat("No exons found in", file ,"and thus no output file was created!\n")
  }
}

