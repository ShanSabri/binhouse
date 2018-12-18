## HELPER FXNS
.current_pkgs <- function() {
    tmp <- installed.packages()
    current_pkgs <- as.vector(tmp[is.na(tmp[, "Priority"]), 1])
    return(current_pkgs)
}

.compare_pkgs <- function(previous_pkgs) {
    missing_pkgs <- setdiff(previous_pkgs, .current_pkgs())
    return(missing_pkgs)
}




## SAVE CURRENT PACKAGES AS RDA
path <- paste0("pkgs_", gsub(" ", "_", R.Version()$version.string))
path <- gsub("\\(|\\)", "", path)
path <- gsub("-", "", path)
assign(path, .current_pkgs())
save(list = path, file = paste0(path, ".rda"))




## ----------------------------------------------------------------
## UPGRADE R
## ----------------------------------------------------------------




## REINSTALL PACKAGES
l <- list.files(pattern = "pkgs")
path <- l[length(l)]
load(file = path)
previous_pkgs <- eval(as.name(gsub(".rda", "", path)))

missing_pkgs <- .compare_pkgs(previous_pkgs)
install.packages(missing_pkgs)
update.packages()

source("https://bioconductor.org/biocLite.R")
biocLite()

missing_pkgs <- .compare_pkgs(previous_pkgs)
for (i in 1:length(missing)) { biocLite(missing_pkgs[i]) }