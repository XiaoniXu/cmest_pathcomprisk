
library(timereg)
library(survival)
library(data.table)
library(mvtnorm)
library(foreach)
library(doParallel)
library(devtools)

# Load the package (which I just fixed)
load_all("cmcrhazard")

# --- Define OLD Method Functions here (to ensure they are correct) ---
getvarnames_old <- function(formula, data = NULL) {
  if (is.character(formula)) return(list(varnames=formula, xvar=formula, yvar=NULL))
  if (is.null(formula)) return(list(varnames=NULL, xvar=NULL, yvar=NULL))
  formula <- formula(formula)
  lyv <- NULL
  lxv <- lvnm <- all.vars(formula[1:2])
  if (length(formula)==3) {
    lyv <- lxv 
    lxv <- all.vars(formula[-2])
    lvnm <- c(lxv, lvnm)
  }
  list(varnames=lvnm, xvar=lxv, yvar=lyv)
}

med_longitudinal_old = function(L=NULL, M, m, Y, treat, control.value, treat.value, data, time_points, peryr=100000) {
  N=dim(data)[1]
  NL=length(L)
  NM=length(M)
  
  MModel = list()
  for (i in 1:NM) MModel[[i]] <- rmvnorm(1, mean = coef(M[[i]]), sigma = vcov(M[[i]]))
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
  
  # Predict M2/L2
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
      if(i==2) pred.data.a.astar.m[, m[1]] <- PredictM_a[,1]
      
      if(NL > 1){
        m1mat.a.m <- model.matrix(~.,data=pred.data.a.m[which(PredictL_a[,i-1]==0),])
        m1mat.astar.m <- model.matrix(~., data = pred.data.astar.m[which(PredictL_astar[,i-1]==0),])
        m1mat.a.astar.m <- model.matrix(~., data = pred.data.a.astar.m[which(PredictL_astar_a[,i-1]==0),])
        PredictM_a[which(PredictL_a[,i-1]==0),i] <- tcrossprod(MModel[[i]], m1mat.a.m)
        PredictM_astar[which(PredictL_astar[,i-1]==0),i] <- tcrossprod(MModel[[i]], m1mat.astar.m)
        PredictM_a_astar[which(PredictL_astar_a[,i-1]==0),i] <- tcrossprod(MModel[[i]], m1mat.a.astar.m)
      } else {
        m1mat.a.m <- model.matrix(~.,data=pred.data.a.m)
        m1mat.astar.m <- model.matrix(~., data = pred.data.astar.m)
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
  pred.data.a.y <- pred.data.astar.y <- pred.data.astar.a.y <- pred.data.astar.a.astar.a.y <- data[,c('idno',getvarnames_old(Y$call)$xvar[-2],m,colnames(data)[grep('time.since.first.exam', colnames(data))])]
  pred.data.a.y[, treat] <- control.value
  pred.data.astar.y[, treat] <- pred.data.astar.a.y[, treat] <- pred.data.astar.a.astar.a.y[, treat] <- treat.value
  pred.data.a.y[, m] <- pred.data.astar.a.y[, m] <- PredictM_a
  pred.data.astar.y[, m] <- PredictM_astar
  pred.data.astar.a.astar.a.y[, m] <- PredictM_a_astar
  pred.data.astar.a.astar.a.y[, m[1]] <- PredictM_a[,1]
  
  vector_time_points <- c()
  for (i in 1:length(time_points)) vector_time_points <- c(vector_time_points, m[i], time_points[i])
  
  f_reshape <- function(df) {
    long <- reshape(df, direction = "long", varying = vector_time_points, sep = "_", times=as.character(seq_along(time_points)), idvar='idno')
    long <- long[order(long$idno),]
    mat <- model.matrix(~., data = long[, match(getvarnames_old(Y$call)$xvar, colnames(long))])[,-1]
    return(mat)
  }
  
  Y_a <- mean(tcrossprod(YModel, f_reshape(pred.data.a.y)))
  Y_astar_a_a <- mean(tcrossprod(YModel, f_reshape(pred.data.astar.a.y)))
  Y_astar_a_astar <- mean(tcrossprod(YModel, f_reshape(pred.data.astar.a.astar.a.y)))
  Y_astar <- mean(tcrossprod(YModel, f_reshape(pred.data.astar.y)))
  
  DE <- (Y_astar_a_a - Y_a) * peryr
  IED <- (Y_astar_a_astar - Y_astar_a_a) * peryr
  IEM <- (Y_astar - Y_astar_a_astar) * peryr
  TE <- (Y_astar - Y_a) * peryr
  return(cbind(DE, IEM, IED, TE))
}

# --- Load data ---
load("simData.mimick.SHS.RDa")
set.seed(123)
n_sub <- 500
df <- simData[sample(nrow(simData), n_sub), ]
df$idno <- 1:nrow(df)
df$M_1 <- df$M1; df$M_2 <- df$M2; df$D_1 <- df$D1; df$D_2 <- df$D2; df$E <- as.numeric(df$E); df$C <- as.numeric(df$C)
df$time.since.first.exam_1 <- 0; df$time.since.first.exam_2 <- df$time_to_event

# --- Fit Original Models ---
mreg <- list(lm(M_1 ~ E + C, data = df), lm(M_2 ~ E + C + M_1, data = df))
D <- list(glm(D_1 ~ E + C + M_1, data = df, family = binomial()), glm(D_2 ~ E + C + M_2, data = df, family = binomial()))
df_process <- tmerge(df, df, id = idno, endpt = event(time_to_event, eventHappened))
df_tv <- reshape(as.data.frame(df), direction = "long", varying = list(M = c("M_1", "M_2"), time = c("time.since.first.exam_1", "time.since.first.exam_2")), 
                 v.names = c("M", "time.since.first.exam"), times = c(1, 2), idvar = "idno")
df_tv <- df_tv[order(df_tv$idno, df_tv$time.since.first.exam), ]
df_process <- tmerge(df_process, df_tv, id = idno, M = tdc(time.since.first.exam, M))
yreg <- aalen(Surv(tstart, tstop, endpt) ~ const(E) + const(M) + const(C), data = df_process, resample.iid = 1, n.sim = 0)
a <- as.numeric(quantile(df$E, 0.25)); astar <- as.numeric(quantile(df$E, 0.75))

# --- Run OLD ---
cat("Running OLD method (re-fitting in loop)...\n")
registerDoParallel(cores = 4)
nboot <- 20
old_boot <- foreach(i=1:nboot, .combine='rbind', .packages=c("mvtnorm", "stats", "survival", "timereg")) %dopar% {
  ind <- sample(1:nrow(df), replace=TRUE)
  data_res <- df[ind,]; data_res$idno <- 1:nrow(data_res)
  mb1 = lm(M_1 ~ E + C, data = data_res); mb2 = lm(M_2 ~ E + C + M_1, data = data_res)
  db1 = glm(D_1 ~ E + C + M_1, data = data_res, family = binomial()); db2 = glm(D_2 ~ E + C + M_2, data = data_res, family = binomial())
  dp_b <- tmerge(data_res, data_res, id = idno, endpt = event(time_to_event, eventHappened))
  dtv_b <- reshape(as.data.frame(data_res), direction = "long", varying = list(M = c("M_1", "M_2"), time = c("time.since.first.exam_1", "time.since.first.exam_2")), 
                   v.names = c("M", "time.since.first.exam"), times = c(1, 2), idvar = "idno")
  dtv_b <- dtv_b[order(dtv_b$idno), ]
  dp_b <- tmerge(dp_b, dtv_b, id = idno, M = tdc(time.since.first.exam, M))
  yb = aalen(Surv(tstart, tstop, endpt) ~ const(E) + const(M) + const(C), data = dp_b, resample.iid = 1, n.sim = 0)
  
  # Map med_longitudinal_old correctly
  pred_res = med_longitudinal_old(L=list(db1, db2), M=list(mb1, mb2), m=c("M_1", "M_2"), Y=yb, treat="E", control.value=a, treat.value=astar, data=data_res, time_points=c("time.since.first.exam_1", "time.since.first.exam_2"))
  return(pred_res)
}
cat("OLD Median Estimates:\n"); print(apply(old_boot, 2, median))

# --- Run NEW (Exact Match mode) ---
cat("\nRunning NEW method (package - EXACT MATCH MODE)...\n")
new_res <- cmest_pathcomprisk(
  D = D, 
  mreg = mreg, 
  mvar = c("M_1", "M_2"), 
  yreg = yreg, 
  avar = "E", 
  a = a, 
  astar = astar, 
  data = df, 
  time_points = c("time.since.first.exam_1", "time.since.first.exam_2"), 
  nboot = nboot, 
  refit = TRUE, 
  yreg_time = "time_to_event", 
  yreg_event = "eventHappened"
)
print(new_res)
