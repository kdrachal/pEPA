

### computes tests of equal predictive accuracy for panels of forecasts 
### as in R. Qu, A. Timmermann and Y. Zhu (2024) 
### Comparing forecasting performance with panel data,
### International Journal of Forecasting 40, 918-941
###
### additionally some tests as in O. Akgun, A. Pirotte, G. Urga and Z. Yang (2024) 
### Equal predictive ability tests based on panel data with applications to OECD and IMF forecasts,
### International Journal of Forecasting 40, 202-228
### are added


##################################################################


### computes loss function L
###
### realized - nxT matrix of realized values
###
### evaluated1 - nxT matrix of forecasts from method 1
###                T - length of time-series
###                n - number of cross-sections
###
### evaluated2 - nxT matrix of forecasts from method 2
###                T - length of time-series
###                n - number of cross-sections
### 
### loss.type - "SE" for squared errors,
###             "AE" for absolute errors,
###             "SPE" for squared proportional error
###                   (useful if errors are heteroskedastic)
###                    see S. J. Taylor, 2005. Asset Price Dynamics, Volatility, and Prediction,
###                        Princeton University Press,
###             "ASE" for absolute scaled error (R. J. Hyndman, A. B. Koehler, 2006,
###                                              Another Look at Measures of Forecast Accuracy,
###                                              International Journal of Forecasting volume 22, 
###                                              679-688, 
###             positive numeric value for loss function of type
###               exp(loss*errors)-1-loss*errors
###               (useful when it is more costly to underpredict y than to overpredict)
###                see U. Triacca,  Comparing Predictive Accuracy of Two Forecasts,
###                    https://www.lem.sssup.it/phd/documents/Lesson19.pdf 


.loss <- function(realized,evaluated1,evaluated2,loss.type)
  {
    e1 <- evaluated1
    e2 <- evaluated2
    for (i in 1:nrow(evaluated1))
      {
        e1[i,] <- as.vector(realized[i,]) - as.vector(evaluated1[i,])
        e2[i,] <- as.vector(realized[i,]) - as.vector(evaluated2[i,])
      }
    
    if (loss.type=="SE") 
      {
        e1 <- e1^2
        e2 <- e2^2
      }
    if (loss.type=="AE") 
      {
        e1 <- abs(e1)
        e2 <- abs(e2)
      }
    if (loss.type=="SPE") 
      {
        for (i in 1:nrow(e1))
          {
            e1[i,] <- (as.vector(e1[i,]) / as.vector(evaluated1[i,]))^2
            e2[i,] <- (as.vector(e2[i,]) / as.vector(evaluated2[i,]))^2
          }
      }
    if (loss.type=="ASE") 
      {
       for (i in 1:nrow(e1))
          {
            temp <- as.vector(realized[i,])
            e1[i,] <- abs(as.vector(e1[i,])) / mean(abs(temp-c(NA,temp[-length(temp)]))[-1])   
            e2[i,] <- abs(as.vector(e2[i,])) / mean(abs(temp-c(NA,temp[-length(temp)]))[-1])   
          }
       e1 <- e1[,-1,drop=FALSE]
       e2 <- e2[,-1,drop=FALSE]
      }
    if (is.numeric(loss.type)) 
      {
        e1 <- exp(loss.type*e1)-1-loss.type*e1
        e2 <- exp(loss.type*e2)-1-loss.type*e2
      }
    
    return(e1-e2)
  }


### test for the pooled average
###
### H0: pooled average loss is equal in expectation for a pair of forecasts from both methods
### H_alt: differences do not average out across the cross-sectional and time-series dimensions
###
### J - maximum lag length


pool_av.test <- function(evaluated1,evaluated2,realized,loss.type="SE",J=NULL)
  {
    nn <- paste(deparse(substitute(evaluated1)),deparse(substitute(evaluated2)),deparse(substitute(realized)),sep=" and ")

    d_L <- .loss(realized,evaluated1,evaluated2,loss.type)

    T <- ncol(d_L)
    n <- nrow(d_L)
    
    d_L_bar <- colMeans(d_L)
    R <- (n^0.5)*d_L_bar
    R_h_bar <- sum(R)/T
    R_bar <- R-R_h_bar

    gammahat <- function(j)
      {
        temp1 <- R_bar
        temp2 <- rep(NA,j)
        temp2 <- c(temp2,temp1)
        temp2 <- temp2[1:T]
        temp <- temp1*temp2
        temp <- temp[(1+j):T]
        temp <- sum(temp)/T
        return(temp)
      }

    if (is.null(J)) { J <- round(T^(1/3)) }
    # U. Triacca,  Comparing Predictive Accuracy of Two Forecasts,
    # https://www.lem.sssup.it/phd/documents/Lesson19.pdf
    
    gdk <- lapply(seq(from=1,to=J,by=1),gammahat)
    gdk <- unlist(gdk)      
    gdk <- c(gdk,gdk)
    
    temp <- seq(from=1,to=J,by=1)
    temp1 <- 1+(temp/J)
    temp2 <- 1-(temp/J)
    temp <- c(temp1,temp2)
    temp <- temp*gdk
    temp <- sum(temp)+gammahat(0)
    temp <- temp^0.5


    J_DM <- sum(d_L)/temp
    J_DM <- J_DM*((n*T)^(-0.5))

    pval <- 2 * min(pnorm(J_DM,lower.tail=FALSE),1 - pnorm(J_DM,lower.tail=FALSE)) 
   
    names(J_DM) <- "statistic"
    names(J) <- "lag length"
    alt <- "Differences between forecasts from method 1 and method 2 do not average out across the cross-sectional and time-series dimensions." 

    ret <- list(J_DM,J,paste(alt),pval,"QTZ test for the pooled average",nn)
    names(ret) <- c("statistic","parameter","alternative","p.value","method","data.name")
    class(ret) <-"htest"
    return(ret)
  }


### test for time clusters
###
### H0: equal predictive accuracy for two methods holds within each of the time clusters
###
### cl - vector of pre-defined blocks of time (starting indices of each block)
###      cl <- c(1, )
###
### NOTICE: The test works if either:
###         a) K >= 2 and significance level <= 0.08326
###         b) 2 <= K <= 14 and significance level <= 0.1
###         c) (K = 2 or K = 3) and significance level <= 0.2


tc.test <- function(evaluated1,evaluated2,realized,loss.type="SE",cl)
  {
    nn <- paste(deparse(substitute(evaluated1)),deparse(substitute(evaluated2)),deparse(substitute(realized)),sep=" and ")

    d_L <- .loss(realized,evaluated1,evaluated2,loss.type)

    T <- ncol(d_L)
    n <- nrow(d_L)
    K <- length(cl)
    
    cl <- c(cl,T+1)
    q <- diff(cl)
    # q <- rep(floor(T/K),K)
    
    Rs <- vector()
    for (i in 1:K)
      {
        d_L_bar <- colMeans(d_L[,((cl[i]):(cl[i+1]-1)),drop=FALSE])
        R <- (n^0.5)*d_L_bar
        R <- sum(R)/q[i]
        Rs[i] <- R
      }
    
    R_bar <- sum(Rs)/K
   
    J_R <- (Rs-R_bar)^2
    J_R <- sum(J_R)/(K-1)
    J_R <- ((K^0.5)*R_bar)/(J_R^0.5)
   
    pval <- 2 * min(pt(q=J_R,df=K-1,lower.tail=FALSE),1 - pt(q=J_R,df=K-1,lower.tail=FALSE))
   
    names(J_R) <- "statistic"
    names(K) <- "number of time clusters"
    alt <- "Equal predictive accuracy for two methods does not hold within each of the time clusters." 

    ret <- list(J_R,K,paste(alt),pval,"QTZ test for for time clusters",nn)
    names(ret) <- c("statistic","parameter","alternative","p.value","method","data.name")
    class(ret) <-"htest"
    return(ret)
  }


### test for cross-sectional clusters
###
### H0: A pair of forecasts have the same expected accuracy among cross-sectional clusters.
###     (Their predictive accuracy could be different across the clusters, but the same among each cluster.)
###
### cl - vector of pre-defined clusters (starting indices of each block)
###      cl <- c(1, )
###
### dc - logical indicating if apply decorrelating clusters
###
### NOTICE: The test works if either:
###         a) K >= 2 and significance level <= 0.08326
###         b) 2 <= K <= 14 and significance level <= 0.1
###         c) (K = 2 or K = 3) and significance level <= 0.2


csc.test <- function(evaluated1,evaluated2,realized,loss.type="SE",cl,dc=FALSE)
  {
    nn <- paste(deparse(substitute(evaluated1)),deparse(substitute(evaluated2)),deparse(substitute(realized)),sep=" and ")

    d_L <- .loss(realized,evaluated1,evaluated2,loss.type)

    T <- ncol(d_L)
    n <- nrow(d_L)
    K <- length(cl)
    
    cl <- c(cl,n+1)
    H <- diff(cl)
    
    if (dc==FALSE)
      {
        Ds <- vector()
        for (i in 1:K)
          {
            D <- d_L[((cl[i]):(cl[i+1]-1)),,drop=FALSE]
            D <- sum(D)
            D <- (H[i]^(-0.5))*(T^(-0.5))*D
            Ds[i] <- D
          }
      }
    else
      {
        As <- matrix(NA,nrow=K,ncol=T)
        for (i in 1:K)
          {
            A <- d_L[((cl[i]):(cl[i+1]-1)),,drop=FALSE]
            A <- colMeans(A)
            As[i,] <- A
          }
        Om <- matrix(0,nrow=K,ncol=K)
        for (i in 1:T)
          {
            A <- As[,i,drop=FALSE]
            A <- A%*%t(A)
            Om <- Om+A
          }
        Om <- Om/T  
        Om <- .sqrtmat(Om)
        Om <- solve(Om)  
        Ds <- as.vector(Om%*%rowSums(As))
      }
   
    D_bar <- sum(Ds)/K
   
    J_D <- (Ds-D_bar)^2
    J_D <- sum(J_D)/(K-1)
    J_D <- ((K^0.5)*D_bar)/(J_D^0.5)
   
    pval <- 2 * min(pt(q=J_D,df=K-1,lower.tail=FALSE),1 - pt(q=J_D,df=K-1,lower.tail=FALSE))
   
    names(J_D) <- "statistic"
    names(K) <- "number of cross-sectional clusters"
    alt <- "A pair of forecasts do not have the same expected accuracy among cross-sectional clusters." 

    ret <- list(J_D,K,paste(alt),pval,"QTZ test for for cross-sectional clusters",nn)
    names(ret) <- c("statistic","parameter","alternative","p.value","method","data.name")
    class(ret) <-"htest"
    return(ret)
  }

.sqrtmat <- function(A)
  {
    if (ncol(A)>1)
      {
        r <- eigen(A)
        V <- r$vectors
        D <- diag(r$values)
        D <- sqrt(D)
    
        out <- V %*% D %*% t(V)
      }
    else
      {
        out <- sqrt(A)
      }
    
    return(out)
  }


### S1 test for the pooled average
###
### H0: pooled average loss is equal in expectation for a pair of forecasts from both methods
### H_alt: differences do not average out across the cross-sectional and time-series dimensions
###
### suitable for situations with cross-sectional independence


pool_av.S1.test <- function(evaluated1,evaluated2,realized,loss.type="SE")
  {
    nn <- paste(deparse(substitute(evaluated1)),deparse(substitute(evaluated2)),deparse(substitute(realized)),sep=" and ")

    d_L <- .loss(realized,evaluated1,evaluated2,loss.type)

    T <- ncol(d_L)
    n <- nrow(d_L)
    
    d_L_bar <- sum(d_L)/(n*T)
   
    d_L_bar_T <- rowMeans(d_L)
    
    d_L_til <- d_L-d_L_bar_T
    
    sig <- 0
    for (i in 1:n)
      {
        for (t in 1:T)
          {
            for (s in 1:T)
              {
                temp1 <- d_L_til[i,t]*d_L_til[i,s]
                temp2 <- (abs(t-s))/(T^(1/3))
               
                if (abs(temp2) <= 1)
                  {
                    temp2 <- 1-abs(temp2)
                  }
                else
                  {
                    temp2 <- 0
                  }
                
                sig <- sig+(temp1*temp2)  
              }
          }
      }
    sig <- sig/(n*T)
    
    S1 <- (((n*T)^0.5)*d_L_bar)/(sig^0.5)

    pval <- 2 * min(pnorm(S1,lower.tail=FALSE),1 - pnorm(S1,lower.tail=FALSE)) 
   
    names(S1) <- "statistic"
    alt <- "Differences between forecasts from method 1 and method 2 do not average out across the cross-sectional and time-series dimensions." 

    ret <- list(S1,paste(alt),pval,"S1 APUY test for the pooled average",nn)
    names(ret) <- c("statistic","alternative","p.value","method","data.name")
    class(ret) <-"htest"
    return(ret)
  }


### S3 test for the pooled average
###
### H0: pooled average loss is equal in expectation for a pair of forecasts from both methods
### H_alt: differences do not average out across the cross-sectional and time-series dimensions
###
### allows for strong cross-sectional dependence


pool_av.S3.test <- function(evaluated1,evaluated2,realized,loss.type="SE")
  {
    nn <- paste(deparse(substitute(evaluated1)),deparse(substitute(evaluated2)),deparse(substitute(realized)),sep=" and ")

    d_L <- .loss(realized,evaluated1,evaluated2,loss.type)

    T <- ncol(d_L)
    n <- nrow(d_L)
    
    d_L_bar <- sum(d_L)/(n*T)
   
    d_L_bar_T <- rowMeans(d_L)
    
    d_L_til <- d_L-d_L_bar_T
    
    sig <- 0
    for (i in 1:n)
      {
        for (j in 1:n)
          {
            for (t in 1:T)
              {
                for (s in 1:T)
                  {
                    temp1 <- d_L_til[i,t]*d_L_til[j,s]
                    temp2 <- (abs(t-s))/(T^(1/3))
                   
                    if (abs(temp2) <= 1)
                      {
                        temp2 <- 1-abs(temp2)
                      }
                    else
                      {
                        temp2 <- 0
                      }
                    
                    sig <- sig+(temp1*temp2)  
                  }
              }
          }
      }
    sig <- sig/(n*T)
    
    S3 <- (((n*T)^0.5)*d_L_bar)/(sig^0.5)

    pval <- 2 * min(pnorm(S3,lower.tail=FALSE),1 - pnorm(S3,lower.tail=FALSE)) 
   
    names(S3) <- "statistic"
    alt <- "Differences between forecasts from method 1 and method 2 do not average out across the cross-sectional and time-series dimensions." 

    ret <- list(S3,paste(alt),pval,"S3 APUY test for the pooled average",nn)
    names(ret) <- c("statistic","alternative","p.value","method","data.name")
    class(ret) <-"htest"
    return(ret)
  }


### C1 test for cross-sectional clusters
###
### H0: A pair of forecasts have the same expected accuracy among cross-sectional clusters.
###     (Their predictive accuracy could be different across the clusters, but the same among each cluster.)
###
### cl - vector of pre-defined clusters (starting indices of each block)
###      cl <- c(1, )
###
### assumes cross-sectional independence


csc.C1.test <- function(evaluated1,evaluated2,realized,loss.type="SE",cl)
  {
    nn <- paste(deparse(substitute(evaluated1)),deparse(substitute(evaluated2)),deparse(substitute(realized)),sep=" and ")

    d_L <- .loss(realized,evaluated1,evaluated2,loss.type)
    
    T <- ncol(d_L)
    n <- nrow(d_L)
    K <- length(cl)
    
    d_L_bar <- sum(d_L)/(n*T)
   
    d_L_bar_T <- rowMeans(d_L)
    
    d_L_til <- d_L-d_L_bar_T

    cl <- c(cl,n+1)
    H <- diff(cl)
    
    id <- c(NA)
    for (i in 1:K)
      {
        id <- c(id,rep(i,H[i]))
      }
    id <- id[-1]
       
    Ds <- vector()
    for (i in 1:K)
      {
        D <- d_L[((cl[i]):(cl[i+1]-1)),,drop=FALSE]
        D <- sum(D)
        D <- (H[i]^(-0.5))*(T^(-0.5))*D
        Ds[i] <- D
      }
 
    I <- diag(1,nrow=K,ncol=K)
 
    sig <- matrix(0,nrow=K,ncol=K)
    for (i in 1:n)
      {
        g <- I[,id[i],drop=FALSE]
        for (t in 1:T)
          {
            for (s in 1:T)
              {
                temp1 <- d_L_til[i,t]*d_L_til[i,s]
                temp2 <- (abs(t-s))/(T^(1/3))
               
                if (abs(temp2) <= 1)
                  {
                    temp2 <- 1-abs(temp2)
                  }
                else
                  {
                    temp2 <- 0
                  }
                
                sig <- sig+(((as.numeric(temp1)*temp2)/H[id[i]])*(g%*%t(g))) 
              }
          }
      }
    sig <- sig/T
 
    C1 <- (t(as.matrix(Ds)))%*%(solve(sig))%*%(as.matrix(Ds))
   
    pval <- 2 * min(pchisq(q=C1,df=K,lower.tail=FALSE),1 - pchisq(q=C1,df=K,lower.tail=FALSE))
   
    names(C1) <- "statistic"
    names(K) <- "number of cross-sectional clusters"
    alt <- "A pair of forecasts do not have the same expected accuracy among cross-sectional clusters." 

    ret <- list(C1,K,paste(alt),pval,"C1 APUY test for for cross-sectional clusters",nn)
    names(ret) <- c("statistic","parameter","alternative","p.value","method","data.name")
    class(ret) <-"htest"
    return(ret)
  }


### C3 test for cross-sectional clusters
###
### H0: A pair of forecasts have the same expected accuracy among cross-sectional clusters.
###     (Their predictive accuracy could be different across the clusters, but the same among each cluster.)
###
### cl - vector of pre-defined clusters (starting indices of each block)
###      cl <- c(1, )
###
### allows for strong cross-sectional dependence


csc.C3.test <- function(evaluated1,evaluated2,realized,loss.type="SE",cl)
  {
    nn <- paste(deparse(substitute(evaluated1)),deparse(substitute(evaluated2)),deparse(substitute(realized)),sep=" and ")

    d_L <- .loss(realized,evaluated1,evaluated2,loss.type)
    
    T <- ncol(d_L)
    n <- nrow(d_L)
    K <- length(cl)
    
    d_L_bar <- sum(d_L)/(n*T)
   
    d_L_bar_T <- rowMeans(d_L)
    
    d_L_til <- d_L-d_L_bar_T

    cl <- c(cl,n+1)
    H <- diff(cl)
    
    id <- c(NA)
    for (i in 1:K)
      {
        id <- c(id,rep(i,H[i]))
      }
    id <- id[-1]
       
    Ds <- vector()
    for (i in 1:K)
      {
        D <- d_L[((cl[i]):(cl[i+1]-1)),,drop=FALSE]
        D <- sum(D)
        D <- (H[i]^(-0.5))*(T^(-0.5))*D
        Ds[i] <- D
      }
 
    I <- diag(1,nrow=K,ncol=K)
 
    sig <- matrix(0,nrow=K,ncol=K)
    for (i in 1:n)
      {
        for (j in 1:n)
          {
            for (t in 1:T)
              {
                for (s in 1:T)
                  {
                    temp1 <- d_L_til[i,t]*d_L_til[i,s]
                    temp2 <- (abs(t-s))/(T^(1/3))
                   
                    if (abs(temp2) <= 1)
                      {
                        temp2 <- 1-abs(temp2)
                      }
                    else
                      {
                        temp2 <- 0
                      }
                    
                    gi <- I[,id[i],drop=FALSE]
                    gj <- I[,id[j],drop=FALSE]
                    sig <- sig+(((as.numeric(temp1)*temp2)/((H[id[i]]*H[id[j]])^0.5))*(gi%*%t(gj))) 
                  }
              }
            }
      }
    sig <- sig/T
 
    C3 <- (t(as.matrix(Ds)))%*%(solve(sig))%*%(as.matrix(Ds))
   
    pval <- 2 * min(pchisq(q=C3,df=K,lower.tail=FALSE),1 - pchisq(q=C3,df=K,lower.tail=FALSE))
   
    names(C3) <- "statistic"
    names(K) <- "number of cross-sectional clusters"
    alt <- "A pair of forecasts do not have the same expected accuracy among cross-sectional clusters." 

    ret <- list(C3,K,paste(alt),pval,"C3 APUY test for for cross-sectional clusters",nn)
    names(ret) <- c("statistic","parameter","alternative","p.value","method","data.name")
    class(ret) <-"htest"
    return(ret)
  }

