# Inspect the mimick SHS data
setwd("c:/Users/berns/Desktop/Research/Valeri/Codes/med_longitudinal_command")
load("simData.mimick.SHS.RDa")

# The object name might be simData or shs_mimick or something else.
# Let's list objects in the environment after loading.
obj_names <- ls()
sink("data_inspection.txt")
cat("Objects loaded:", obj_names, "\n\n")

for (obj in obj_names) {
    cat("--- Object:", obj, "---\n")
    cat("Class:", class(get(obj)), "\n")
    if (is.data.frame(get(obj))) {
        cat("Dimensions:", dim(get(obj)), "\n")
        cat("Column names:\n")
        print(colnames(get(obj)))
        cat("\nHead:\n")
        print(head(get(obj)))
    } else {
        print(get(obj))
    }
    cat("\n")
}
sink()
