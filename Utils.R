
#' gffReads reads in the annotation portion of a GFF3 file and returns a data 
#' frame containing the contents.
#'
#' @param gffFile File location of GFF3 file
#' @param nrows the number of rows to read in the GFF3 file (this is important to use
#' when skipping the fasta portion of the file - which will prevent this from running
#' properly)
#'
#' @return a data frame containing the sequence name, source, feature
#' start, end, score, strand, frame, and attributes columns.
#'
#' @examples  df <- gffRead(file, nrows = fasta_line - comment_num - 1)
gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, " ... ", sep="")
  gff <- read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows)
  names(gff) <- c("seqname", "source", "feature", "start", "end",
                 "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows")
  return(gff)
}


#' getAttributeField separates attribute field from complete attributes string
#'
#' @param x string version of the attributes column 
#' @param field field to be separated
#' @param attrsep (default ";") character separating fields in attributes section
#'
#' @return a string containing the desired attribute separated from the others
#'
#' @examples df$Name <- getAttributeField(df$attributes, "Name")
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}


