
# Functions from 250422_SHS_gfor_single_exp.R

getvarnames <- function(formula, data = NULL)
{
  if (is.character(formula))
    return(list(varnames=formula, xvar=formula, yvar=NULL))
  if (is.null(formula)) return(list(varnames=NULL, xvar=NULL, yvar=NULL))
  
  formula <- formula(formula)
  lyv <- NULL
  lxv <- lvnm <- all.vars(formula[1:2])
  if (length(formula)==3) {
    lyv <- lxv 
    lxv <- all.vars(formula[-2])
    if ("." %in% lxv) {
      if (length(data)==0)
        stop("!getvarnames! '.' in formula and no 'data'")
      lform <- formula(terms(formula, data=data))
      lxv <- all.vars(lform[-2])
    }
    lvnm <- c(lxv, lvnm)
  }
  list(varnames=lvnm, xvar=lxv, yvar=lyv)
}

med_longitudinal=function(L=NULL, M, m, Y, treat='logcocr', control.value=0, treat.value=1, data, time_points, peryr=100000){
  
  N=dim(data)[1]
  NL=length(L)
  NM=length(M)
  
  MModel = list()
  for (i in 1:NM){
    MModel[[i]] <- mvtnorm::rmvnorm(1, mean = coef(M[[i]]), sigma = vcov(M[[i]]))
  }
  
  YModel = mvtnorm::rmvnorm(1, mean = Y$gamma, sigma = Y$robvar.gamma)
  
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
  
  DE <- mean(PredictY_DEIEM - PredictY_TEDE_2)*peryr
  IEM <- mean(PredictY_IEDTE_1 - PredictY_IEMIED)*peryr
  IED <- mean(PredictY_IEMIED - PredictY_DEIEM)*peryr
  TE <- mean(PredictY_IEDTE_1 - PredictY_TEDE_2)*peryr
  effects <- cbind(DE, IEM, IED, TE)
  return(effects)
}
