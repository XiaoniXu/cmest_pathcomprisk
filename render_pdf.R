# Script to sanitize walkthrough.md and render to PDF
library(rmarkdown)

# Paths
input_file <- "C:/Users/berns/.gemini/antigravity/brain/fa0ff944-1aa3-4fdc-9a90-2095922cad80/walkthrough.md"
sanitized_file <- "walkthrough_sanitized.md"
output_pdf <- "c:/Users/berns/Desktop/Research/Valeri/Codes/med_longitudinal_command/walkthrough.pdf"

# Read content
content <- readLines(input_file, warn = FALSE)

# Sanitize Emojis
content <- gsub("\u2705", "[X]", content) # Check mark
content <- gsub("\u26A0\uFE0F", "[!]", content) # Warning
content <- gsub("\u2139\uFE0F", "[i]", content) # Info
content <- gsub("\ud83d\udca1", "[TIP]", content) # Light bulb
content <- gsub("\u2757", "[!]", content) # Exclamation

# Sanitize file links (remove file:///c:/path/to/ and just keep the basename or relative path)
# e.g. [test_cmest_shs.R](file:///c:/Users/berns/Desktop/Research/Valeri/Codes/med_longitudinal_command/test_cmest_shs.R)
# to [test_cmest_shs.R](test_cmest_shs.R)
content <- gsub("file:///c:/Users/berns/Desktop/Research/Valeri/Codes/med_longitudinal_command/", "", content)

# Write sanitized content
writeLines(content, sanitized_file)

# Render to PDF
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
rmarkdown::render(sanitized_file, output_format = "pdf_document", output_file = output_pdf)
