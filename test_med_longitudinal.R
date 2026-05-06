# Refactored test script using cmcrhazard package
library(devtools)
library(timereg)
library(tidyr)
library(dplyr)
library(foreach)
library(data.table)
library(survival)
library(mvtnorm)
library(doParallel)

# Set up parallel processing
registerDoParallel(cores = 4)

# Set working directory
setwd("c:/Users/berns/Desktop/Research/Valeri/Codes/med_longitudinal_command")

# Load the package
devtools::load_all("cmcrhazard")

# 1. Load simulated data
load("simData.mimick.SHS.RDa")

# 2. Data Preparation
# Using full dataset (no sampling)
simData <- data.table::as.data.table(simData)


# Prepare variables
simData[, M_1 := M1]
simData[, M_2 := M2]
simData[, D_1 := D1]
simData[, D_2 := D2]

# Unique ID and format
simData[, idno := .I]
df_main <- as.data.frame(simData)

# Counting-process helpers (required for outcome model)
df_main$time_to_event[df_main$time_to_event <= 0] <- 0.0001
df_main$time.since.first.exam_1 <- 0
df_main$time.since.first.exam_2 <- df_main$time_to_event

# 3. Model Fitting
cat("Fitting models...\n")

# Mediator models (mreg)
mreg1 <- lm(M_1 ~ E + C, data = df_main)
mreg2 <- lm(M_2 ~ E + M_1 + C, data = df_main)
mreg <- list(mreg1, mreg2)

# Dropout/Survival models (D)
D1_mod <- glm(D_1 ~ E + M_1 + C, data = df_main, family = binomial())
D2_mod <- glm(D_2 ~ E + M_2 + C, data = df_main, family = binomial())
dreg <- list(D1_mod, D2_mod)

# Outcome model (yreg - Aalen)
# Must use the same naming convention for M as in the package (reshaped to 'M')
cat("Preparing counting-process data for Aalen model...\n")
df_process <- tmerge(df_main, df_main, id = idno, endpt = event(time_to_event, eventHappened))

# Reshape mediators for tdc()
df_tv <- reshape(df_main,
    direction = "long",
    varying = list(
        M = c("M_1", "M_2"),
        time = c("time.since.first.exam_1", "time.since.first.exam_2")
    ),
    v.names = c("M", "time.since.first.exam"),
    times = c(1, 2), idvar = "idno"
)
df_tv <- df_tv[order(df_tv$idno, df_tv$time.since.first.exam), ]
df_process <- tmerge(df_process, df_tv, id = idno, M = tdc(time.since.first.exam, M))

# df_tv <- reshape(data, direction = "long", varying = c("sbp_1", "time.since.first.exam_1", "sbp_2", "time.since.first.exam_2", "sbp_3", "time.since.first.exam_3"), sep = "_", times=c('1','2','3'), idvar='id_unique')

# Fit Aalen model
yreg <- aalen(Surv(tstart, tstop, endpt) ~ const(E) + const(M) + const(C),
    data = df_process, resample.iid = 1, n.sim = 0
)

# 4. Mediation Analysis using the Package
cat("Running mediation analysis with the cmcrhazard package...\n")
cat(paste0("Bootstrap samples: 1000 (n = ", nrow(df_main), ")\n\n"))

# Median values for treatment comparison
a <- as.numeric(quantile(df_main$E, 0.25, na.rm = TRUE))
astar <- as.numeric(quantile(df_main$E, 0.75, na.rm = TRUE))

start_time <- Sys.time()

results <- cmest_pathcomprisk(
    dreg = dreg,
    mreg = mreg,
    mvar = c("M_1", "M_2"),
    yreg = yreg,
    avar = "E",
    a = a,
    astar = astar,
    data = df_main,
    time_points = c("time.since.first.exam_1", "time.since.first.exam_2"),
    peryr = 100000,
    nboot = 1000,
    refit = TRUE,                # Enable Exact Match (Non-parametric Bootstrap)
    yreg_time = "time_to_event", # Time-to-event variable
    yreg_event = "eventHappened" # Event indicator variable
)

length_time <- Sys.time() - start_time
cat("Time taken:", round(length_time, 2), "units\n")
print(results)

# Save results
saveRDS(results, "med_longitudinal_package_results.rds")
cat("\nResults saved to: med_longitudinal_package_results.rds\n")
