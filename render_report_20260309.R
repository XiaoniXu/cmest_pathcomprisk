
library(rmarkdown)

# Specify the path to Pandoc found on this system
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")

# Render the RMarkdown to PDF
rmarkdown::render("reports/20260309_Progress_Report_XX.Rmd", output_format = "pdf_document")
