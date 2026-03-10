
library(rmarkdown)

# Paths
input_file <- "c:/Users/berns/Desktop/Research/Valeri/Codes/med_longitudinal_command/20260223_Progress_Report_XX.Rmd"
output_pdf <- "c:/Users/berns/Desktop/Research/Valeri/Codes/med_longitudinal_command/20260223_Progress_Report_XX.pdf"

# Set Pandoc path (from render_pdf.R)
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")

# Render
rmarkdown::render(input_file, output_format = "pdf_document", output_file = output_pdf)
