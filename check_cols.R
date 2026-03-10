load("simData.Y10D10.Population.RData")
if (exists("simData")) {
    print(colnames(simData))
    print(nrow(simData))
} else {
    print("simData not found")
}
