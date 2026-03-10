# Debug script for med_longitudinal function

library(SimDesign)
library(timereg)
library(tidyr)
library(dplyr)
library(foreach)
library(data.table)
library(survival)
library(mvtnorm)

# Run sequentially for debugging
foreach::registerDoSEQ()

# Set working directory
setwd("c:/Users/berns/Desktop/Research/Valeri/Codes/med_longitudinal_command")

# Load simulated data
load("simData.Y10D10.Population.RData")

# Subset to very small sample for debugging
set.seed(123)
simData <- simData[sample(1:nrow(simData), size = 100, replace = FALSE), ]

# Prepare data
simData[, M_1 := M1]
simData[, M_2 := M2]

# Counting-process helpers
simData[, `:=`(
    tstart = 0, tstop = time_to_event,
    time.since.first.exam_1 = 0,
    time.since.first.exam_2 = time_to_event
)]

# Unique ID
simData[, id_unique := .I]
simData <- as.data.frame(simData)

cat("Data prepared. Sample size:", nrow(simData), "\n")
cat("Running debug analysis with 1 sample...\n\n")

getvarnames <- function(formula, data = NULL) {
    if (is.character(formula)) {
        return(list(varnames = formula, xvar = formula, yvar = NULL))
    }
    if (is.null(formula)) {
        return(list(varnames = NULL, xvar = NULL, yvar = NULL))
    }
    formula <- formula(formula)
    lyv <- NULL
    lxv <- lvnm <- all.vars(formula[1:2])
    if (length(formula) == 3) {
        lyv <- lxv
        lxv <- all.vars(formula[-2])
        if ("." %in% lxv) {
            if (length(data) == 0) {
                stop("!getvarnames! '.' in formula and no 'data'")
            }
            lform <- formula(terms(formula, data = data))
            lxv <- all.vars(lform[-2])
        }
        lvnm <- c(lxv, lvnm)
    }
    list(varnames = lvnm, xvar = lxv, yvar = lyv)
}

# Main med_longitudinal function with debug prints
med_longitudinal <- function(D = NULL, M, m, Y, treat = "logcocr", control.value = a, treat.value = a_star, data, time_points, peryr = 100000) {
    cat("Entering med_longitudinal...\n")
    N <- dim(data)[1]
    ND <- length(D)
    NM <- length(M)
    cat("N:", N, "ND:", ND, "NM:", NM, "\n")

    MModel <- list()
    for (i in 1:NM) {
        MModel[[i]] <- rmvnorm(1, mean = coef(M[[i]]), sigma = vcov(M[[i]]))
    }

    YModel <- rmvnorm(1, mean = Y$gamma, sigma = Y$robvar.gamma)

    PredictD_a <- PredictD_astar <- PredictD_astar_a <- matrix(NA, nrow = N, ncol = ND)
    PredictM_a <- PredictM_astar <- PredictM_a_astar <- matrix(NA, nrow = N, ncol = NM)

    # Predict M1
    cat("Predicting M1...\n")
    pred.data.astar.m1 <- pred.data.a.m1 <- model.frame(M[[1]])
    pred.data.astar.m1[, treat] <- treat.value
    pred.data.a.m1[, treat] <- control.value

    m1mat.astar <- model.matrix(terms(M[[1]]), data = pred.data.astar.m1)
    m1mat.a <- model.matrix(terms(M[[1]]), data = pred.data.a.m1)

    PredictM_astar[, 1] <- tcrossprod(MModel[[1]], m1mat.astar)
    PredictM_a[, 1] <- tcrossprod(MModel[[1]], m1mat.a)

    # Predict D1
    if (ND > 0) {
        cat("Predicting D1...\n")
        pred.data.astar.d1 <- pred.data.a.d1 <- pred.data.astar.a.d1 <- model.frame(D[[1]])
        pred.data.astar.d1[, treat] <- pred.data.astar.a.d1[, treat] <- treat.value
        pred.data.a.d1[, treat] <- control.value
        pred.data.astar.d1[, m[1]] <- PredictM_astar[, 1]
        pred.data.a.d1[, m[1]] <- pred.data.astar.a.d1[, m[1]] <- PredictM_a[, 1]

        PredictD_astar_a[, 1] <- rbinom(N, size = 1, prob = predict(D[[1]], pred.data.astar.a.d1, type = "response"))
        PredictD_a[, 1] <- rbinom(N, size = 1, prob = predict(D[[1]], pred.data.a.d1, type = "response"))
        PredictD_astar[, 1] <- rbinom(N, size = 1, prob = predict(D[[1]], pred.data.astar.d1, type = "response"))
    }

    # Predict Li/Di
    if (NM > 1) {
        for (i in 2:NM) {
            cat("Predicting M", i, "...\n")
            pred.data.a.m <- pred.data.astar.m <- pred.data.a.astar.m <- as.data.frame(matrix(nrow = N, ncol = (dim(model.frame(M[[i]]))[2] - 1)))
            colnames(pred.data.a.m) <- colnames(pred.data.astar.m) <- colnames(pred.data.a.astar.m) <- attr(terms(M[[i]]), "term.labels")
            names <- colnames(pred.data.a.m)[which(colnames(pred.data.a.m) %in% attr(terms(M[[1]]), "term.labels"))]

            pred.data.a.m[, names] <- pred.data.a.astar.m[, names] <- pred.data.astar.m[, names] <- model.frame(M[[1]])[, names]
            pred.data.a.m[, treat] <- pred.data.a.astar.m[, treat] <- control.value
            pred.data.astar.m[, treat] <- treat.value

            pred.data.a.m[, m[i - 1]] <- PredictM_a[, i - 1]
            pred.data.astar.m[, m[i - 1]] <- PredictM_astar[, i - 1]
            pred.data.a.astar.m[, m[i - 1]] <- PredictM_a_astar[, i - 1]
            if (i == 2) {
                pred.data.a.astar.m[, m[1]] <- PredictM_a[, 1]
            }

            if (ND > 1) {
                cat("Filtering by previous D...\n")
                m1mat.a.m <- model.matrix(~., data = pred.data.a.m[which(PredictD_a[, i - 1] == 0), ])
                m1mat.astar.m <- model.matrix(~., data = pred.data.astar.m[which(PredictD_astar[, i - 1] == 0), ])
                m1mat.a.astar.m <- model.matrix(~., data = pred.data.a.astar.m[which(PredictD_astar_a[, i - 1] == 0), ])

                PredictM_a[which(PredictD_a[, i - 1] == 0), i] <- tcrossprod(MModel[[i]], m1mat.a.m)
                PredictM_astar[which(PredictD_astar[, i - 1] == 0), i] <- tcrossprod(MModel[[i]], m1mat.astar.m)
                PredictM_a_astar[which(PredictD_astar_a[, i - 1] == 0), i] <- tcrossprod(MModel[[i]], m1mat.a.astar.m)
            } else {
                m1mat.a.m <- model.matrix(~., data = pred.data.a.m)
                m1mat.astar.m <- model.matrix(~., data = pred.data.astar.m)
                m1mat.a.astar.m <- model.matrix(~., data = pred.data.a.astar.m)

                PredictM_a[, i] <- tcrossprod(MModel[[i]], m1mat.a.m)
                PredictM_astar[, i] <- tcrossprod(MModel[[i]], m1mat.astar.m)
                PredictM_a_astar[, i] <- tcrossprod(MModel[[i]], m1mat.a.astar.m)
            }

            if (ND > 1 & i <= ND) {
                cat("Predicting D", i, "...\n")
                pred.data.a.d <- pred.data.astar.d <- pred.data.astar.a.d <- as.data.frame(matrix(nrow = N, ncol = (dim(model.frame(D[[i]]))[2] - 1)))
                colnames(pred.data.a.d) <- colnames(pred.data.astar.d) <- colnames(pred.data.astar.a.d) <- attr(terms(D[[i]]), "term.labels")

                # ... skipping some assignments for brevity, focusing on crash location
                names <- colnames(pred.data.a.d)[which(colnames(pred.data.a.d) %in% attr(terms(D[[1]]), "term.labels"))]

                pred.data.a.d[, names] <- pred.data.astar.a.d[, names] <- pred.data.astar.d[, names] <- model.frame(D[[1]])[, names]
                pred.data.a.d[, treat] <- control.value
                pred.data.astar.d[, treat] <- pred.data.astar.a.d[, treat] <- treat.value

                pred.data.a.d[, m[i]] <- PredictM_a[, i]
                pred.data.astar.d[, m[i]] <- PredictM_astar[, i]
                pred.data.astar.a.d[, m[i]] <- PredictM_a_astar[, i]

                cat("Taking rbinom for D", i, "...\n")
                PredictD_a[which(PredictD_a[, i - 1] == 0), i] <- rbinom(length(which(PredictD_a[, i - 1] == 0)), size = 1, prob = predict(D[[i]], pred.data.a.d[which(PredictD_a[, i - 1] == 0), ], type = "response"))
                PredictD_astar[which(PredictD_astar[, i - 1] == 0), i] <- rbinom(length(which(PredictD_astar[, i - 1] == 0)), size = 1, prob = predict(D[[i]], pred.data.astar.d[which(PredictD_astar[, i - 1] == 0), ], type = "response"))
                PredictD_astar_a[which(PredictD_astar_a[, i - 1] == 0), i] <- rbinom(length(which(PredictD_astar_a[, i - 1] == 0)), size = 1, prob = predict(D[[i]], pred.data.astar.a.d[which(PredictD_astar_a[, i - 1] == 0), ], type = "response"))
            }
        }
    }

    cat("Success so far...\n")
    return(matrix(1, 4, 1))
}


i <- 1
cat("Loop iteration 1\n")
ind <- sample(1:nrow(simData), replace = TRUE)
data <- simData[ind, ]
data$id_unique <- 1:nrow(data)
data$idno <- if ("id" %in% names(data)) data$id else data$id_unique

if (!("M" %in% names(data))) data$M <- data$M_1

M1 <- lm(M_1 ~ E + C, data = data)
L1 <- glm(D1 ~ E + M_1 + C, data = data, family = binomial())
M2 <- lm(M_2 ~ E + M_1 + C, data = data)
L2 <- glm(D2 ~ E + M_2 + C, data = data, family = binomial())

df <- tmerge(data, data, id = id_unique, endpt = event(time_to_event, eventHappened))

df_tv <- reshape(as.data.frame(data),
    direction = "long",
    varying = list(
        M = c("M_1", "M_2"),
        time = c("time.since.first.exam_1", "time.since.first.exam_2")
    ),
    v.names = c("M", "time.since.first.exam"),
    times = c(1, 2), idvar = "id_unique"
)
df_tv <- df_tv[order(df_tv$id_unique, df_tv$time), ]
df <- tmerge(df, df_tv, id = id_unique, M = tdc(time.since.first.exam, M))

Y <- aalen(Surv(tstart, tstop, endpt) ~ const(E) + const(M) + const(C),
    data = df, resample.iid = 1, n.sim = 0
)

treat <- "E"
D <- list(D1 = L1, D2 = L2)
Mlist <- list(M1 = M1, M2 = M2)
mvec <- c("M_1", "M_2")
tpts <- c("1", "2")

a <- as.numeric(quantile(data$E, 0.25, na.rm = TRUE))
a_star <- as.numeric(quantile(data$E, 0.75, na.rm = TRUE))

cat("Calling med_longitudinal...\n")
out <- med_longitudinal(
    D = D, M = Mlist, m = mvec, Y = Y, treat = treat,
    control.value = a, treat.value = a_star,
    data = data, time_points = tpts, peryr = 100000
)
print(out)
