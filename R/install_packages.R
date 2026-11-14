# Installs what the analysis needs. INLA is not on CRAN and comes from its own repository.
pkgs <- c("sf", "Matrix", "dplyr", "tidyr", "ggplot2", "knitr")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing)

if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages("INLA",
                   repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
                   dep = TRUE)
}

library(INLA)
cat("INLA", as.character(packageVersion("INLA")), "is ready\n")
