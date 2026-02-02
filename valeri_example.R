# Example script to run cmest_pathcomprisk with simulated data
# Updated to reflect new function name and arguments

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

# 3. Load Simulated Data
load("simData.mimick.SHS.RDa")

# 4. Data Pre-processing
set.seed(2)
simData <- simData[sample(1:nrow(simData), size = 200, replace = FALSE), ] 

df <- as.data.frame(simData)
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

df$eventHappened <- as.integer(df$eventHappened > 0)
if ("D1" %in% names(df)) df$D1 <- as.integer(df$D1 > 0)
if ("D2" %in% names(df)) df$D2 <- as.integer(df$D2 > 0)

keep <- complete.cases(df[, c("E","C","M_1","M_2")])
df2  <- df[keep, , drop = FALSE]

# 5. Define Models

# mreg (Mediator models) - formerly M
mreg <- list(
  lm(M_1 ~ E + C,           data = df2),
  lm(M_2 ~ E + M_1 + C,     data = df2)
)

# D (Competing Risk models) - formerly L
D <- list(
  glm(D1 ~ E + M_1 + C, data = df2, family = binomial()),
  glm(D2 ~ E + M_2 + C, data = df2, family = binomial())
)

# yreg (Outcome model) - formally Y
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

# Arguments updated:
# L -> D
# M -> mreg
# m -> mvar
# Y -> yreg
# treat -> avar
# control.value -> a
# treat.value -> astar
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
