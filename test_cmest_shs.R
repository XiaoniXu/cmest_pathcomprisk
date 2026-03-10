# Analysis script for simData.mimick.SHS.RDa using cmcrhazard package

# 1. Environment Setup
library(timereg)
library(survival)
library(data.table)
library(mvtnorm)
library(foreach)
library(doParallel)

# Find R executable path (hardcoded from previous discovery)
R_PATH <- "C:/Program Files/R/R-4.5.2/bin/x64/Rscript.exe"

# Install the package if not already installed or to ensure latest version
# setwd("c:/Users/berns/Desktop/Research/Valeri/Codes/med_longitudinal_command")
# install.packages("./cmcrhazard", repos = NULL, type = "source")
# devtools::install("./cmcrhazard")
devtools::load_all("cmcrhazard")

# Set up parallel processing
registerDoParallel(cores = 4)

# 2. Data Loading & Preparation
load("simData.mimick.SHS.RDa")
# Subsample for efficiency (5000 rows as per plan)
set.seed(123)
simData <- simData[sample(1:nrow(simData), size = 5000, replace = FALSE), ]

# Prepare variables
simData$idno <- simData$id
simData$M_1 <- simData$M1
simData$M_2 <- simData$M2
simData$D_1 <- simData$D1
simData$D_2 <- simData$D2
# Create time columns matching the package's expectations
simData$time.since.first.exam_1 <- 0
simData$time.since.first.exam_2 <- simData$time_to_event

# 3. Model Fitting
cat("Fitting models...\n")

# Mediator models
# M1 depends on E and C
mreg1 <- lm(M_1 ~ E + C, data = simData)
# M2 depends on E, C, and previous M (M1)
# Note: The function expects the mediator names to match the wide format
mreg2 <- lm(M_2 ~ E + C + M_1, data = simData)
mreg <- list(mreg1, mreg2)

# Covariate/Dropout models
# D1 depends on E, C, and M1
D1_mod <- glm(D_1 ~ E + C + M_1, data = simData, family = binomial())
# D2 depends on E, C, and M2
D2_mod <- glm(D_2 ~ E + C + M_2, data = simData, family = binomial())
D <- list(D1_mod, D2_mod)

# Outcome model (Aalen additive hazards)
# We need to prepare counting-process data for the Aalen model fit
# In the long format used by cmest_pathcomprisk, the mediator is named 'M'
# because it reshapes M_1, M_2 with sep='_'.

cat("Preparing counting-process data for Aalen model...\n")
# Corrected tmerge call:
# Initial merge to set up endpt
df_process <- tmerge(simData, simData, id = idno, endpt = event(time_to_event, eventHappened))

# Reshape mediators into long format for tdc()
df_tv <- reshape(as.data.frame(simData),
    direction = "long",
    varying = list(
        M = c("M_1", "M_2"),
        time = c("time.since.first.exam_1", "time.since.first.exam_2")
    ),
    v.names = c("M", "time.since.first.exam"),
    times = c(1, 2), idvar = "idno"
)
df_tv <- df_tv[order(df_tv$idno, df_tv$time.since.first.exam), ]

# Apply tdc to get time-varying M
df_process <- tmerge(df_process, df_tv, id = idno, M = tdc(time.since.first.exam, M))

# Fit the Aalen model
# Important: Use column names that will match the long-format created by cmest_pathcomprisk
yreg <- aalen(Surv(tstart, tstop, endpt) ~ const(E) + const(M) + const(C),
    data = df_process, resample.iid = 1, n.sim = 0
)

# 4. Mediation Analysis
cat("Running mediation analysis with cmest_pathcomprisk...\n")
cat("Using 100 bootstrap samples...\n")

# Median values for treatment comparison (a vs astar)
a <- as.numeric(quantile(simData$E, 0.25))
astar <- as.numeric(quantile(simData$E, 0.75))

results <- cmest_pathcomprisk(
    D = D,
    mreg = mreg,
    mvar = c("M_1", "M_2"),
    yreg = yreg,
    avar = "E",
    a = a,
    astar = astar,
    data = simData,
    time_points = c("time.since.first.exam_1", "time.since.first.exam_2"),
    peryr = 100000,
    nboot = 100,
    refit = TRUE,                # Enable Exact Match
    yreg_time = "time_to_event", # Time variable
    yreg_event = "eventHappened" # Event variable
)

# 5. Output Results
cat("\nResults for mimick SHS Data:\n")
print(results)

# Save results
saveRDS(results, "shs_mimick_results.rds")
cat("\nResults saved to shs_mimick_results.rds\n")
