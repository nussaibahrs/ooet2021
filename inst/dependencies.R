# Attachments ----
to_install <- c("caret", "chronosphere", "data.table", "divDyn", "doParallel", "dplyr", "foreach", "forecast", "ggplot2", "ggpubr", "ggthemes", "janitor", "magrittr", "mgcv", "nlme", "patchwork", "purrr", "raster", "RColorBrewer", "reshape2", "scales", "sp", "tidyverse", "viridis")
  for (i in to_install) {
    message(paste("looking for ", i))
    if (!requireNamespace(i)) {
      message(paste("     installing", i))
      install.packages(i)
    }
  }
