library(SimDesign)
library(timereg)
library(tidyr)
library(dplyr)
library(mice)
library(splitstackshape)
library(foreach)
library(mvtnorm)
# Only if you want to run parallel models
#library(doParallel)

#num_cores <- parallel::detectCores() - 2
#cl <- makeCluster(num_cores)
#registerDoParallel(cl)
# Ensure the cluster is stopped when the function/script exits
#on.exit(stopCluster(cl))




load("simData.mimick.SHS.RDa")

# Mediation main function
med_longitudinal=function(L=NULL, M, m, Y, treat='logcocr', control.value=a, treat.value=a_star, data, time_points, peryr=100000){
  
  N=dim(data)[1]
  NL=length(L)
  NM=length(M)
  
  MModel = list()
  for (i in 1:NM){
    MModel[[i]] <- rmvnorm(1, mean = coef(M[[i]]), sigma = vcov(M[[i]]))
  }
  
  YModel = rmvnorm(1, mean = Y$gamma, sigma = Y$robvar.gamma)
  
  PredictL_a <- PredictL_astar <- PredictL_astar_a <- matrix(NA,nrow=N,ncol=NL)
  PredictM_a <- PredictM_astar <- PredictM_a_astar <- matrix(NA,nrow=N,ncol=NM)
  
  # Predict M1
  pred.data.astar.m1 <- pred.data.a.m1 <- model.frame(M[[1]])
  pred.data.astar.m1[, treat] <- treat.value
  pred.data.a.m1[, treat] <- control.value
  
  m1mat.astar <- model.matrix(terms(M[[1]]), data = pred.data.astar.m1)
  m1mat.a <- model.matrix(terms(M[[1]]), data = pred.data.a.m1)
  
  PredictM_astar[,1] <- tcrossprod(MModel[[1]], m1mat.astar)
  PredictM_a[,1] <- tcrossprod(MModel[[1]], m1mat.a)
  
  # Predict L1
  if(NL > 0){
    pred.data.astar.l1 <- pred.data.a.l1 <- pred.data.astar.a.l1 <- model.frame(L[[1]])
    pred.data.astar.l1[, treat] <- pred.data.astar.a.l1[, treat] <- treat.value
    pred.data.a.l1[, treat] <- control.value
    pred.data.astar.l1[, m[1]] <- PredictM_astar[,1]
    pred.data.a.l1[, m[1]] <- pred.data.astar.a.l1[, m[1]] <- PredictM_a[,1]
    
    PredictL_astar_a[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.astar.a.l1, type='response'))
    PredictL_a[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.a.l1, type='response'))
    PredictL_astar[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.astar.l1, type='response'))
  }
  
  # Predict Li (only if more than one time point)
  if (NM > 1){
    for (i in 2:NM){
      pred.data.a.m <- pred.data.astar.m <- pred.data.a.astar.m <- as.data.frame(matrix(nrow=N, ncol=(dim(model.frame(M[[i]]))[2]-1)))
      colnames(pred.data.a.m) <- colnames(pred.data.astar.m) <- colnames(pred.data.a.astar.m) <- attr(terms(M[[i]]),"term.labels")
      names <- colnames(pred.data.a.m)[which(colnames(pred.data.a.m) %in% attr(terms(M[[1]]),"term.labels"))]
      
      pred.data.a.m[, names] <- pred.data.a.astar.m[, names] <- pred.data.astar.m[, names] <- model.frame(M[[1]])[,names]
      pred.data.a.m[, treat] <- pred.data.a.astar.m[, treat] <- control.value
      pred.data.astar.m[, treat] <- treat.value
      
      pred.data.a.m[, m[i-1]] <- PredictM_a[,i-1]
      pred.data.astar.m[, m[i-1]] <- PredictM_astar[,i-1]
      pred.data.a.astar.m[, m[i-1]] <- PredictM_a_astar[,i-1]
      if(i==2){
        pred.data.a.astar.m[, m[1]] <- PredictM_a[,1]
      }
      
      if(NL > 1){
        m1mat.a.m <- model.matrix(~.,data=pred.data.a.m[which(PredictL_a[,i-1]==0),])
        m1mat.astar.m <- model.matrix(~., data = pred.data.astar.m[which(PredictL_astar[,i-1]==0),])
        m1mat.a.astar.m <- model.matrix(~., data = pred.data.a.astar.m[which(PredictL_astar_a[,i-1]==0),])
        
        PredictM_a[which(PredictL_a[,i-1]==0),i] <- tcrossprod(MModel[[i]], m1mat.a.m)
        PredictM_astar[which(PredictL_astar[,i-1]==0),i] <- tcrossprod(MModel[[i]], m1mat.astar.m)
        PredictM_a_astar[which(PredictL_astar_a[,i-1]==0),i] <- tcrossprod(MModel[[i]], m1mat.a.astar.m)
      } else{
        m1mat.a.m <- model.matrix(~.,data=pred.data.a.m)
        m1mat.astar.m <- model.matrix(~., data = pred.data.astar.m)
        m1mat.a.astar.m <- model.matrix(~., data = pred.data.a.astar.m)
        
        PredictM_a[,i] <- tcrossprod(MModel[[i]], m1mat.a.m)
        PredictM_astar[,i] <- tcrossprod(MModel[[i]], m1mat.astar.m)
        PredictM_a_astar[,i] <- tcrossprod(MModel[[i]], m1mat.a.astar.m)
      }
      
      if(NL > 1 & i<=NL){
        pred.data.a.l <- pred.data.astar.l <- pred.data.astar.a.l <- as.data.frame(matrix(nrow=N, ncol=(dim(model.frame(L[[i]]))[2]-1)))
        colnames(pred.data.a.l) <- colnames(pred.data.astar.l) <- colnames(pred.data.astar.a.l) <- attr(terms(L[[i]]),"term.labels")
        names <- colnames(pred.data.a.l)[which(colnames(pred.data.a.l) %in% attr(terms(L[[1]]),"term.labels"))]
        
        pred.data.a.l[, names] <- pred.data.astar.a.l[, names] <- pred.data.astar.l[, names] <- model.frame(L[[1]])[,names]
        pred.data.a.l[, treat] <- control.value
        pred.data.astar.l[, treat] <- pred.data.astar.a.l[, treat] <- treat.value
        
        pred.data.a.l[, m[i]] <- PredictM_a[,i]
        pred.data.astar.l[, m[i]] <- PredictM_astar[,i]
        pred.data.astar.a.l[, m[i]] <- PredictM_a_astar[,i]
        
        PredictL_a[which(PredictL_a[,i-1]==0),i] <- rbinom(length(which(PredictL_a[,i-1]==0)), size=1, prob=predict(L[[i]], pred.data.a.l[which(PredictL_a[,i-1]==0),], type='response'))
        PredictL_astar[which(PredictL_astar[,i-1]==0),i] <- rbinom(length(which(PredictL_astar[,i-1]==0)), size=1, prob=predict(L[[i]], pred.data.astar.l[which(PredictL_astar[,i-1]==0),], type='response'))
        PredictL_astar_a[which(PredictL_astar_a[,i-1]==0),i] <- rbinom(length(which(PredictL_astar_a[,i-1]==0)), size=1, prob=predict(L[[i]], pred.data.astar.a.l[which(PredictL_astar_a[,i-1]==0),], type='response'))
      }
    }
  }
  
  
  # Predict Y
  # Data augmentation method for person-time database
  # PredictY_DEIEM: a*, D1_a, M1_aD1a, D2_aD1aM1aD1a, M2_aD1aM1aD2a
  # PredictY_TEDE_2: a, D1_a, M1_aD1a, D2_aD1a M1aD1a, M2_aD1aM1aD2a
  # PredictY_IEMIED: a*, D1_a, M1_a*D1a, D2_aD1a M1a*D1a, M2_a*D1aM1a*D2a
  # PredictY_IEDTE_1: a*, D1_a*, M1_a*D1a*, D2_a*D1a*M1a*D1a*, M2_a*D1a*M1a*D2a*
  
  pred.data.a.y <- pred.data.astar.y <- pred.data.astar.a.y <- pred.data.astar.a.astar.a.y <- data[,c('id_unique',getvarnames(Y$call)$xvar[-2],m,colnames(data)[grep('time.since.first.exam', colnames(data))])]
  pred.data.a.y[, treat] <- control.value
  pred.data.astar.y[, treat] <- pred.data.astar.a.y[, treat] <- pred.data.astar.a.astar.a.y[, treat] <- treat.value
  pred.data.a.y[, m] <- pred.data.astar.a.y[, m] <- PredictM_a
  pred.data.astar.y[, m] <- PredictM_astar
  pred.data.astar.a.astar.a.y[, m] <- PredictM_a_astar
  pred.data.astar.a.astar.a.y[, m[1]] <- PredictM_a[,1]
  
  
  ########################
  
  # Data augmentation method for the counterfactuals
  vector_time_points <- c()
  for (i in 1:length(time_points)){
    vector_time_points <- c(vector_time_points, m[i], time_points[i])
  }
  
  # pred.data.a.y
  df_tv <- reshape(pred.data.a.y, direction = "long", varying = vector_time_points,
                   sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df_pred.data.a.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
  df_pred.data.a.y <- model.matrix(~.,data=df_pred.data.a.y)[,-1]
  
  # pred.data.astar.y
  df_tv <- reshape(pred.data.astar.y, direction = "long", varying = vector_time_points,
                   sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df_pred.data.astar.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
  df_pred.data.astar.y <- model.matrix(~.,data=df_pred.data.astar.y)[,-1]
  
  # pred.data.astar.a.astar.a.y
  df_tv <- reshape(pred.data.astar.a.astar.a.y, direction = "long", varying = vector_time_points,
                   sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df_pred.data.astar.a.astar.a.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
  df_pred.data.astar.a.astar.a.y <- model.matrix(~.,data=df_pred.data.astar.a.astar.a.y)[,-1]
  
  # pred.data.astar.a.y
  df_tv <- reshape(pred.data.astar.a.y, direction = "long", varying = vector_time_points,
                   sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df_pred.data.astar.a.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
  df_pred.data.astar.a.y <- model.matrix(~.,data=df_pred.data.astar.a.y)[,-1]
  
  #######################
  
  PredictY_DEIEM <- mean(tcrossprod(YModel, df_pred.data.astar.a.y))
  PredictY_TEDE_2 <- mean(tcrossprod(YModel, df_pred.data.a.y))
  PredictY_IEMIED <- mean(tcrossprod(YModel, df_pred.data.astar.a.astar.a.y))
  PredictY_IEDTE_1 <- mean(tcrossprod(YModel, df_pred.data.astar.y))
  
  DE <- mean(PredictY_DEIEM - PredictY_TEDE_2)*100000
  IEM <- mean(PredictY_IEDTE_1 - PredictY_IEMIED)*100000
  IED <- mean(PredictY_IEMIED - PredictY_DEIEM)*100000
  TE <- mean(PredictY_IEDTE_1 - PredictY_TEDE_2)*100000
  effects <- cbind(DE, IEM, IED, TE)
  return(effects)
}


set.seed(2)  
# subset the data into 10000 instead of 2000. 2025/8/26
simData <- simData[sample(1:nrow(simData), size = 10000, replace = FALSE), ]



# ---- Start: data prep for simData ----
df <- as.data.frame(simData)

# 1) Rename & add required columns (wide format expected)
df$idno <- df$id
df$M_1  <- df$M1
df$M_2  <- df$M2

# time index helpers required by your function’s reshape

# simulated time: 4 and 8 years after
df$time.since.first.exam_1 <- 4
df$time.since.first.exam_2 <- 8

# coerce types if needed
to_num <- c("E","C","M_1","M_2","time_to_event")
for (cc in intersect(to_num, names(df))) df[[cc]] <- suppressWarnings(as.numeric(as.character(df[[cc]])))
df$eventHappened <- as.integer(df$eventHappened > 0)
if ("D1" %in% names(df)) df$D1 <- as.integer(as.character(df$D1))
if ("D2" %in% names(df)) df$D2 <- as.integer(as.character(df$D2))

# keep only rows complete for BOTH mediators and predictors they use
keep <- complete.cases(df[, c("E","C","M_1","M_2")])
df2  <- df[keep, , drop = FALSE]   # this will be ~1795 rows 

cat("Rows kept:", nrow(df2), "\n")  #  should print 1795




df2 <- within(df, {
  idno <- id
  id_unique <- id
  M_1  <- M1
  M_2  <- M2
  time.since.first.exam_1 <- 1
  time.since.first.exam_2 <- 2
})
df2 <- df2[complete.cases(df2[, c("M_1","M_2")]), ]

m <- c("M_1","M_2")
time_points <- c("time.since.first.exam_1","time.since.first.exam_2")

control.value <- as.numeric(quantile(df2$E, 0.25, na.rm = TRUE))
treat.value   <- as.numeric(quantile(df2$E, 0.75, na.rm = TRUE))
treat <- "E"

M <- list(
  lm(M_1 ~ E + C,           data = df2),
  lm(M_2 ~ E + M_1 + C,     data = df2)
)

L <- list(
  glm(D1 ~ E + M_1 + C, data = df2, family = binomial()),
  glm(D2 ~ E + M_2 + C, data = df2, family = binomial())
)

base <- df2[, c("idno","E","C","M_1","M_2","time_to_event","eventHappened")]

# Build a long(person-time) data with exactly two rows per subject
longA <- transform(base,
                   tstart = 0, tstop = pmin(1, time_to_event),
                   endpt  = as.integer(eventHappened == 1 & time_to_event <= 1),
                   M      = M_1)
longB <- transform(base,
                   tstart = 1, tstop = time_to_event,
                   endpt  = as.integer(eventHappened == 1 & time_to_event > 1),
                   M      = M_2)

## Keep valid intervals only
longA <- longA[longA$tstop > longA$tstart, ]
longB <- longB[longB$tstop > longB$tstart, ]
df_long <- rbind(longA, longB)

## Fit additive hazards model
Y <- aalen(Surv(tstart, tstop, endpt) ~ const(E) + const(M) + const(C),
           data = df_long, n.sim = 0, resample.iid = 1)

stopifnot(all(c("idno", treat, m,
                "time.since.first.exam_1","time.since.first.exam_2")
              %in% names(df2)))

res <- med_longitudinal(
  L = L,             # or NULL if you don't model dropout
  M = M,
  m = m,
  Y = Y,
  treat = treat,
  control.value = control.value,
  treat.value   = treat.value,
  data = df2,
  time_points = time_points,
  peryr = 100000
)
print(res)