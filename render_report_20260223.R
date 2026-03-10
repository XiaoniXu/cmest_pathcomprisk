
library(rmarkdown)
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
rmarkdown::render("reports/20260223_Progress_Report_XX.Rmd", output_format = "pdf_document")
