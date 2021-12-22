
prediction <- function(data_example = long, surv = surv){
  
  long_example <- data_example[,c("ID", "time", "CEA", "CA199", "CA125")]#7
  
  long_example$CEA_cut <- ifelse(long_example$CEA > 50, 50, long_example$CEA)
  long_example$CA199_cut <- ifelse(long_example$CA199 > 370, 370, long_example$CA199)
  long_example$CA125_cut <- ifelse(long_example$CA125 > 350, 350, long_example$CA125)
  
  long_example$ln_CEA <- log(long_example$CEA_cut)
  long_example$ln_CA199_cut <- log(long_example$CA199_cut)
  long_example$ln_CA125_cut <- log(long_example$CA125_cut)
  
  
  long_example <- long_example[,c("ID", "time", "ln_CEA", "ln_CA199_cut", "ln_CA125_cut")]#7
  
  surv <- unique(surv)#1
  
  long_zero <- surv[, c("ID", "JX_CEA", "JX_CA199", "JX_CA125")]
  long_zero$time <- 0
  long_zero$ln_CEA <- log(min(long_zero$JX_CEA, 50))
  long_zero$ln_CA199_cut <- log(min(long_zero$JX_CA199),370)
  long_zero$ln_CA125_cut <- log(min(long_zero$JX_CA125),350)
  
  long_zero <- long_zero[, c("ID", "time", "ln_CEA", "ln_CA199_cut", "ln_CA125_cut")]#1
  
  long_example <- rbind(long_zero, long_example)#7
  
  time <- long_example$time
  time <- time[2:length(time)]
  
  #载入模型
  load("./list_save.RData")#51
  
  #构建模型
  # RSF Model
  formula_OS <- list_save[["formula_OS"]]
  train.surv.temp <- list_save[["train.surv.temp"]]
  rsf.fit.OS <- rfsrc(as.formula(formula_OS), data = train.surv.temp, ntree = 1000, seed = 1234)
  
  
  pred_marker <- list()
  pred_prob <- list()
  
  #载入rsf模型
  for (j in 1:length(time)) {
    
    t <- time[j]
    tmp.data = long_example[long_example$time <= t,]
    
    ###univariate FPCA in test
    pace_test <- list(pace_CEA = NULL, pace_CA199= NULL, pace_CA125 = NULL)
    
    
    ID_long_test <- tmp.data$ID
    tVec_test <- tmp.data$time
    yVec_CEA_test <- tmp.data$ln_CEA 
    yVec_CA199_test <- tmp.data$ln_CA199_cut 
    yVec_CA125_test <- tmp.data$ln_CA125_cut 
    yVec_three_test <- list(yVec_CEA_test, yVec_CA199_test, yVec_CA125_test)
    
    Xi.test =  predict_test = predgrid = NULL
    
    uPACE <- list_save[["uPACE"]]
    L <- list_save[["L"]]
    
    for (i in 1:3) {
      pace_test[[i]] <- MakeFPCAInputs(
        IDs = ID_long_test,
        tVec = tVec_test,
        yVec = yVec_three_test[[i]],
        na.rm = FALSE,
        sort = FALSE,
        deduplicate = TRUE
      )
      tmp.ufpca <- predict(uPACE[[i]], pace_test[[i]]$Ly, pace_test[[i]]$Lt, K=L[i])
      predgrid <- tmp.ufpca$predGrid
      Xi.test = cbind(Xi.test, tmp.ufpca$scores) # FPC scores
      predict_test[[i]] <- tmp.ufpca
    }
    

    # calculate marker level
    # estimate MFPC scores for test subjects
    Cms <- list_save[["Cms"]]
    rho.test = mfpca.score(Xi.test, Cms) 
    
    # predict longitudinal trajectories 
    meanFun.train <- list_save[["meanFun.train"]]
    psi <- list_save[["psi"]]
    long.pred = mfpca.pred(rho.test, meanFun.train, psi)
    
    workgrid <- uPACE[[1]][["workGrid"]]
    
    CEA_pred <- long.pred[,,1]
    CEA <- as.data.frame(CEA_pred)
    CEA <- exp(CEA)
    
    
    CA199_pred <- long.pred[,,2]
    CA199 <- as.data.frame(CA199_pred)
    CA199 <- exp(CA199)
    
    CA125_pred <- long.pred[,,3]
    CA125 <- as.data.frame(CA125_pred)
    CA125 <- exp(CA125)
    
    marker <- data.frame(t =  rep(t, 51), time = workgrid, 
                         CEA = CEA[,1], CA199 = CA199[,1], CA125 = CA125[,1]) 
    
    pred_marker[[j]] <- marker
    
    
    #FPC scores
    
    # estimate MFPC scores for test subjects 
    rho.test = mfpca.score(Xi.test, Cms)  
    colnames(rho.test) = paste0("rho.", (1: ncol(rho.test)))
    test.surv.temp = cbind(surv[,2:15], rho.test)
    
    deltaT <- seq(0, 60, 0.1)

    DP.prob.rsf.OS = NULL 
    ith=0
    # prediction for different time windowes 
    for(dt in deltaT){
      ith=ith+1
      DP.prob = cond.prob.pec.rsf(rsf.fit.OS, test.surv.temp, t, (t+dt))
      DP.prob.rsf.OS <- c(DP.prob.rsf.OS, DP.prob)
    }
    prob <- data.frame(t =  rep(t, 601), dt = deltaT, surv_prob = DP.prob.rsf.OS) 
    prob$tp <- prob$t + prob$dt  
    
    pred_prob[[j]] <- prob
  } 
  list_output <- list(long_example = long_example, 
                      pred_marker = pred_marker,
                      pred_prob = pred_prob)
  return(list_output)
}



plot <- function(long_example = long_example, pred_marker = pred_marker, pred_prob = pred_prob){
  time <- long_example$time
  time <- time[2:length(time)]
  
  saveGIF({
    for (j in 1:length(time)) {
      
      t <- time[j]
      prob <- pred_prob[[j]]
      marker <- pred_marker[[j]]
      
      true_marker_t <- subset(long_example, long_example$time <= t)
      marker <- subset(marker, marker$time <= t)
      
      dynamic_CEA <- ggplot(prob, aes(x = prob$tp, y = rescale(prob[,3], c(0,60))))+
        geom_smooth(colour = "red", size = 1.2)+
        geom_line(aes(x = marker$time, y = marker[,3]), data = marker, colour = "blue", size = 1)+
        geom_point(aes(x = true_marker_t$time, y = exp(true_marker_t$ln_CEA)), data = true_marker_t, size = 2, pch = 8)+
        scale_y_continuous(limits = c(0, 64), breaks = seq(0, 64, 5), sec.axis = sec_axis(~./60, breaks = c(0,0.2, 0.4,0.6,0.8,1.0), name = "Survival  Probability"))+
        scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 6), sec.axis=sec_axis(~., breaks = NULL))+
        geom_vline(aes(xintercept= t),linetype=5,col="black")+
        theme_classic()+
        labs( x = NULL,
              y = "CEA, ng/ml")+
        theme(axis.title= element_text(size=13, color="black", family="Arial",face = "bold"))+
        theme(axis.text.x= element_text(size=12, color="black", family="Arial",face = "bold"))+
        theme(axis.text.y= element_text(size=12, color="black", family="Arial",face = "bold"))+
        theme(panel.grid = element_blank()) + 
        theme(legend.position = "none")
      print(dynamic_CEA)
    }    
  },ani.width = 700, ani.height=500)
}



plot_end <- function(long_example = long_example, pred_marker = pred_marker, pred_prob = pred_prob){
  time <- long_example$time
  time <- time[2:length(time)]
  
  j <- length(time)    
  t <- time[j]
  prob <- pred_prob[[j]]
  marker <- pred_marker[[j]]
  
  true_marker_t <- subset(long_example, long_example$time <= t)
  marker <- subset(marker, marker$time <= t)
  
  pred_plot <- ggplot(prob, aes(x = prob$tp, y = rescale(prob[,3], c(0,60))))+
    geom_smooth(colour = "red", size = 1.2)+
    geom_line(aes(x = marker$time, y = marker[,3]), data = marker, colour = "blue", size = 1)+
    geom_point(aes(x = true_marker_t$time, y = exp(true_marker_t$ln_CEA)), data = true_marker_t, size = 2, pch = 8)+
    scale_y_continuous(limits = c(0, 64), breaks = seq(0, 64, 5), sec.axis = sec_axis(~./60, breaks = c(0,0.2, 0.4,0.6,0.8,1.0), name = "Survival  Probability"))+
    scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 6), sec.axis=sec_axis(~., breaks = NULL))+
    geom_vline(aes(xintercept= t),linetype=5,col="black")+
    theme_classic()+
    labs( x = NULL,
          y = "CEA, ng/ml")+
    theme(axis.title= element_text(size=13, color="black", family="Arial",face = "bold"))+
    theme(axis.text.x= element_text(size=12, color="black", family="Arial",face = "bold"))+
    theme(axis.text.y= element_text(size=12, color="black", family="Arial",face = "bold"))+
    theme(panel.grid = element_blank()) + 
    theme(legend.position = "none")
  
  return(pred_plot)
}




# multivariate FPCA based on results from uPACE
mFPCA = function(Xi, phi, p , L, I=I){
  
  # eigenanalysis on matrix M
  M = t(Xi) %*% Xi/(I-1)
  eigen.M = eigen(M)
  values = eigen.M$values
  pve = cumsum(values)/sum(values)
  Cms = eigen.M$vectors
  index = unlist(lapply(1:length(L), function(x) rep(x, L[x])))
  
  # MFPCA score
  rho = mfpca.score(Xi, Cms)
  
  # MFPCA eigenfunction
  psis = NULL
  for(j in 1:p){
    psi = NULL
    for(m in 1:dim(Cms)[2]){
      psi = cbind(psi, phi[[j]] %*% Cms[which(index==j),m])
    }
    psis[[j]] = psi
  }
  
  out = list(eigenvalue = values, Cms = Cms, pve = pve, index=index, rho = rho, psis=psis)
  
  return(out)
}

# mfpc score calculation
mfpca.score = function(predXi, Cms){
  rho = matrix(NA, nrow = nrow(predXi), ncol=dim(Cms)[2])
  for(i in 1:nrow(predXi)){
    for(m in 1:dim(Cms)[2]){
      rho[i,m] = predXi[i,] %*% Cms[,m]
    }
  }
  return(rho)
}


# mfpc trajectories prediction
mfpca.pred = function(score, meanf, psi, n.rho=NULL){
  p = length(psi)
  n = nrow(score)
  
  if(is.null(n.rho)){
    n.rho = ncol(score)
  }
  
  pred = array(NA, c(n, length(meanf[[1]]), p))
  for(m in 1:p){
    pred[,,m] = matrix(meanf[[m]], nrow=n, ncol =length(meanf[[m]]), byrow = T ) + score[,1:n.rho]%*%t(psi[[m]][, 1:n.rho])
  }
  
  out = pred
  return(out)
}

# risk predict using predictSurvProb instead of survfit to handle RSF.(REQUIRES PEC)
cond.prob.pec.rsf = function(model, newdata, Tstart, Tpred){
  risk.Tstart = as.numeric(predictSurvProb(model, newdata=newdata, times=Tstart))
  risk.Tpred = as.numeric(predictSurvProb(model, newdata=newdata, times=Tpred))
  return(risk.Tpred/risk.Tstart)
}