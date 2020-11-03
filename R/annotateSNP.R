#' Annotate significant SNPs
#' 
#' Annotation of significant SNPs with gene info.
#' 
#' @export
annotateSNP <- function(GWAS,
                        gff) {
  signSnp <- do.call(rbind, GWAS$signSnp)
  
  gffDat <- data.table::fread(gff, skip = 1, na.strings = c("###", "."), 
                              fill = TRUE, sep = "\t")
  colnames(gffDat) <- c("seqid", "source", "type", "start", "end", "score",
                      "strand", "phase", "attributes")
  gffDat <- gffDat[!is.na(gffDat[["seqid"]]), ]
  gffDat <- gffDat[gffDat[["type"]] == "gene"]
  
  gffDat[["chr"]] <- suppressWarnings(as.numeric(gsub(pattern = "Chr", 
                                                      replacement = "",
                                                      x = gffDat[["seqid"]])))
  gffDat <- gffDat[!is.na(gffDat[["chr"]]), ]
  
  data.table::setalloccol(gffDat, 2048)
  annoDat <- gffDat[signSnp, on = c("chr", "start <= pos", "end >= pos"), 
                    nomatch = NA]
  
  subAttribs <- strsplit(annoDat[["attributes"]], split = ";")
  annoDat[["ID"]] <- sapply(X = subAttribs, FUN = function(subAttrib) {
    if (!is.na(subAttrib[1])) {
      items <- strsplit(subAttrib, split = "=")
      items <- setNames(sapply(items, `[[`, 2), sapply(items, `[[`, 1))
      items[["ID"]]
    } else {
      NA
    }
  })
  
  annoDat2 <- split(annoDat[["ID"]], annoDat[["snp"]])
  
  return(annoDat2)
}

