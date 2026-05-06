library(timereg)
library(survival)
library(data.table)
library(mvtnorm)
library(foreach)
library(doParallel)

# Set up parallel processing
registerDoParallel(cores = 4)

# Load the old functions
source("reports/old_functions.R")

# 1. Load data
load("simData.mimick.SHS.RDa")

# Subsample for efficiency
set.seed(123)
simData <- as.data.frame(simData[sample(1:nrow(simData), size = 5000, replace = FALSE), ])

# 2. Add required exact formatting
# 2a. id_unique must be set
simData$id_unique <- 1:nrow(simData)
simData$M_1 <- simData$M1
simData$M_2 <- simData$M2
simData$D_1 <- simData$D1
simData$D_2 <- simData$D2
simData$E <- as.numeric(simData$E)
simData$C <- as.numeric(simData$C)
simData$time.since.first.exam_1 <- 0
simData$time.since.first.exam_2 <- simData$time_to_event

# 3. Model fitting
cat("Fitting baseline models...\n")
mreg <- list()
mreg[[1]] <- lm(M_1 ~ E + C, data = simData)
mreg[[2]] <- lm(M_2 ~ E + C + M_1, data = simData)

D <- list()
D[[1]] <- glm(D_1 ~ E + C + M_1, data = simData, family = binomial())
D[[2]] <- glm(D_2 ~ E + C + M_2, data = simData, family = binomial())

cat("Preparing counting process data for Aalen...\n")
df_process <- tmerge(simData, simData, id = id_unique, endpt = event(time_to_event, eventHappened))

# reshape to long for time-varying M
df_tv <- reshape(as.data.frame(simData), direction = "long", 
                 varying = list(M = c("M_1", "M_2"), time = c("time.since.first.exam_1", "time.since.first.exam_2")), 
                 v.names = c("M", "time.since.first.exam"), times = c(1, 2), idvar = "id_unique")
df_tv <- df_tv[order(df_tv$id_unique, df_tv$time.since.first.exam), ]

df_process <- tmerge(df_process, df_tv, id = id_unique, M = tdc(time.since.first.exam, M))

# Fit Aalen outcome model
# CRITICAL: const(M) must be the 2nd covariate because of hardcoded xvar[-2]
yreg <- aalen(Surv(tstart, tstop, endpt) ~ const(E) + const(M) + const(C), 
              data = df_process, resample.iid = 1, n.sim = 0)

a <- as.numeric(quantile(simData$E, 0.25))
astar <- as.numeric(quantile(simData$E, 0.75))

# 4. Monte Carlo simulation using the raw med_longitudinal function
cat("Running Monte Carlo with the old med_longitudinal function...\n")

nboot <- 10
start_time <- Sys.time()

# Base R loop
results_mc <- lapply(1:nboot, function(b) {
  res <- med_longitudinal(
    L = D, 
    M = mreg, 
    m = c("M_1", "M_2"), 
    Y = yreg, 
    treat = "E", 
    control.value = a, 
    treat.value = astar, 
    data = simData, 
    time_points = c("time.since.first.exam_1", "time.since.first.exam_2"),
    peryr = 100000
  )
  return(res)
})

results_mc <- do.call(rbind, results_mc)

# 5. Summarize the results
# Calculate mean and 95% CIs
calc_summary <- function(vec) {
  est <- mean(vec, na.rm = TRUE)
  ci_lower <- quantile(vec, 0.025, na.rm = TRUE)
  ci_upper <- quantile(vec, 0.975, na.rm = TRUE)
  paste0(round(est, 2), " (", round(ci_lower, 2), ", ", round(ci_upper, 2), ")")
}

cat("\nSummary of MC run (", nboot, " iterations):\n")
cat("Direct Effect (DE):", calc_summary(results_mc[, 1]), "\n")
cat("Indirect Effect Mediator (IEM):", calc_summary(results_mc[, 2]), "\n")
cat("Indirect Effect Dropout (IED):", calc_summary(results_mc[, 3]), "\n")
cat("Total Effect (TE):", calc_summary(results_mc[, 4]), "\n")

end_time <- Sys.time()
cat(sprintf("Time taken: %.2f mins\n", as.numeric(difftime(end_time, start_time, units = "mins"))))
