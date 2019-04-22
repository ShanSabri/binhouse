### Generate the Normalizer Values for each Sample
get_normalizer <- function(mat, genes.norm.factor = c("ACTB", "CLTC", "RPLP0")) {

  message("Generating the geometric mean of housekeeper genes")

  if (!length(genes.norm.factor) > 1) {
    stop("Must have more than 1 housekeeper gene")
  }

  if (!all(genes.norm.factor %in% rownames(mat))) {
    missing.ind <- !genes.norm.factor %in% rownames(mat)
    stop(paste0("Matrix is missing the genes: ", genes.norm.factor[missing.ind]))
  }

  normalizer.log <- log(mat[genes.norm.factor, ])
  normalizer.log.mean <- apply(normalizer.log, 2, mean)
  normalizer.log.mean.exp <- exp(normalizer.log.mean)
  normalizer.log.mean.exp
}

### Normalize the Expression Matrix 
hk_normalize <- function(mat, normalizer) {

  message("Normalizing the expression matrix")
  mat.norm <- sweep(mat, 2, normalizer, "/")
  mat.norm.scaled <- mat.norm * 1000

  # Ensure that any values < 1 do not become negative
  mat.norm.scaled[mat.norm.scaled < 1] <- 1

  message("Log2 transforming")
  mat.norm.scaled.log <- log(mat.norm.scaled, 2)
  mat.norm.scaled.log
} 

### Example
# data <- .....
normalizer <- get_normalizer(data)
norm <- hk_normalize(data, normalizer)
