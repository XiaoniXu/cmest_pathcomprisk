# Target inspection of simData
setwd("c:/Users/berns/Desktop/Research/Valeri/Codes/med_longitudinal_command")
load("simData.mimick.SHS.RDa")

sink("simData_details.txt")
if (exists("simData")) {
    cat("Class of simData:", class(simData), "\n")
    cat("Dimensions:", dim(simData), "\n")
    cat("Column Names:\n")
    print(colnames(simData))
    cat("\nSummary:\n")
    print(summary(simData))
    cat("\nHead:\n")
    print(head(simData))
} else {
    cat("simData object not found in the file.\n")
}
sink()
