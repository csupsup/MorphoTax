## MorphoTax
This R package contains a collection of functions designed to provide a transparent and reproducible pipeline for analyzing morphological data in taxonomic studies.

### Installation
```
## Check and install required packages
req.pkgs <- c("devtools", "ggplot2", "tibble", "dplyr", "car", "mda",
                  "psych", "magrittr", "caret, "RColorBrewer")

for (pkg in req.pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

## Install "MorphoTax"
install_github("csupsup/MorphoTax")
```
### Using the package
In progress ...
