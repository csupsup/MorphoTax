## MorphoTax: A set of functions compiled to analyze morphological data for taxonomic studies
### Description

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
