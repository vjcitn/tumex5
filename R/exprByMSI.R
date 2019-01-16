#' bind microsatellite instability scores via MSIsensor to a TCGA expression dataset
#' @param tumcode character(1) e.g., "COAD"
#' @param genesym character(1) gene symbol
#' @param alias character(1) an alternative name for gene
#' @note This function works specifically with the *rnagn (RESTful)
#' SummarizedExperiment instance in package tumex5.  The assay
#' data are in hsdshdflab.hdfgroup.org and internet access
#' is require for this function to succeed.
#' @export
 exprByMSI = function(tumcode, genesym, alias) {
  if (missing(alias)) alias=genesym
  ob = paste0(tumcode, "rnagn")
  ex = get(ob)
  ex = BiocOncoTK::bindMSI(ex)
  data.frame(
   patient_barcode=colnames(ex),
   acronym=tumcode,
   symbol = genesym,
   alias = alias,
   log2ex=log2(as.numeric(SummarizedExperiment::assay(ex[genesym,]))+1),
   msicode = ifelse(ex$msiTest >= 4, ">=4", "<4"))
 }

