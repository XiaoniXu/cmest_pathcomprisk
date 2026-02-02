# Example script to run cmest_pathcomprisk with SELF-GENERATED data
# This script simulates data from scratch instead of loading an RData file.
# Uses base R for simulation to avoid dependency issues.

# 1. Load required libraries
library(devtools)
library(timereg)
library(tidyr)
library(dplyr)
library(mice)
library(splitstackshape)
library(foreach)
library(doParallel)

# 2. Load the local package
pkg_path <- "cmcrhazard" 
cat("Loading package from", pkg_path, "\n")
devtools::load_all(pkg_path)

# ==========================================
# 3. Simulate Data (Base R)
# ==========================================
cat("Simulating data without simstudy...\n")
set.seed(456)
n <- 200

# Params
M_0 = 100; M_1_E = 1.2; M_1_C = 3.9
M_2_0 = 46; M_2_E = 0.9; M_2_C = -1.3; M_2_M_1 = 0.5
D_1_0 = -4.8; D_1_E = .28; D_1_C = .62; D_1_M_1 = -.001
Y_1_0 = -7.4; Y_1_E = 0.14; Y_1_C = 0.48; Y_1_M_1 = 0.02

# Generate covariates and mediators
E <- rnorm(n, mean=2, sd=sqrt(0.7))
C <- rbinom(n, size=1, prob=0.4)
M1 <- rnorm(n, mean = M_0 + M_1_E * E + M_1_C * C, sd=1)
M2 <- rnorm(n, mean = M_2_0 + M_2_M_1 * M1 + M_2_E * E + M_2_C * C, sd=1)

# Generate Survival Times (Exponential approximation for simplicity)
# Hazard = exp(LinearPredictor)
# Time = -log(U) / Hazard
lambda_Y <- exp(Y_1_0 + Y_1_E * E + Y_1_C * C + Y_1_M_1 * M1)
lambda_D <- exp(D_1_0 + D_1_E * E + D_1_C * C + D_1_M_1 * M1)

Y_time <- -log(runif(n)) / lambda_Y
D_time <- -log(runif(n)) / lambda_D

simData <- data.frame(E, C, M1, M2, Y=Y_time, D=D_time)
simData$id <- 1:n

# Define indicators based on cutoffs
t_2 = 4
t_3 = 8 

simData$D1 = ifelse((simData$D <= t_2) & (simData$D <= simData$Y), 1, 0)
simData$D2 = ifelse((simData$D > t_2) & (simData$D <= t_3) & (simData$D <= simData$Y), 1, 0)

simData$time_to_event = pmin(simData$Y, simData$D)
simData$eventHappened = ifelse(simData$Y < simData$D, 1, 0)

# Censor at t_3
simData$eventHappened[simData$time_to_event > t_3] = 0
simData$time_to_event[simData$time_to_event > t_3] = t_3


# ==========================================
# 4. Data Pre-processing (Same as valeri_example.R)
# ==========================================

df <- simData
df$idno <- df$id
df$M_1  <- df$M1
df$M_2  <- df$M2

# Create time-indexing columns required by the method
df$time.since.first.exam_1 <- 4
df$time.since.first.exam_2 <- 8

# Ensure variables are numeric/integer as expected by models
to_num <- c("E","C","M_1","M_2","time_to_event", "D1", "D2")
for (cc in intersect(to_num, names(df))) {
  if(is.factor(df[[cc]])) {
    df[[cc]] <- as.numeric(as.character(df[[cc]]))
  } else {
    df[[cc]] <- as.numeric(as.character(df[[cc]]))
  }
}

keep <- complete.cases(df[, c("E","C","M_1","M_2")])
df2  <- df[keep, , drop = FALSE]

# 5. Define Models

# mreg (Mediator models)
mreg <- list(
  lm(M_1 ~ E + C,           data = df2),
  lm(M_2 ~ E + M_1 + C,     data = df2)
)

# D (Competing Risk models)
D <- list(
  glm(D1 ~ E + M_1 + C, data = df2, family = binomial()),
  glm(D2 ~ E + M_2 + C, data = df2, family = binomial())
)

# yreg (Outcome model)
base <- df2[, c("idno","E","C","M_1","M_2","time_to_event","eventHappened")]
longA <- transform(base,
                   tstart = 0, tstop = pmin(1, time_to_event),
                   endpt  = as.integer(eventHappened == 1 & time_to_event <= 1),
                   M      = M_1)
longB <- transform(base,
                   tstart = 1, tstop = time_to_event,
                   endpt  = as.integer(eventHappened == 1 & time_to_event > 1),
                   M      = M_2)
longA <- longA[longA$tstop > longA$tstart, ]
longB <- longB[longB$tstop > longB$tstart, ]
df_long <- rbind(longA, longB)

yreg <- aalen(Surv(tstart, tstop, endpt) ~ const(E) + const(M) + const(C),
              data = df_long, n.sim = 0, resample.iid = 1)


# 6. Run the G-Formula Function
registerDoParallel(cores = 4) 

mvar <- c("M_1","M_2")
time_points <- c("time.since.first.exam_1", "time.since.first.exam_2") 

cat("Running cmest_pathcomprisk...\n")
results <- cmest_pathcomprisk(
  D = D,
  mreg = mreg,
  mvar = mvar,
  yreg = yreg,
  avar = "E",
  a = as.numeric(quantile(df2$E, 0.25, na.rm = TRUE)),
  astar = as.numeric(quantile(df2$E, 0.75, na.rm = TRUE)),
  data = df2,
  time_points = time_points,
  peryr = 100000,
  nboot = 20
)

print(results)
