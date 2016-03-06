#' Monte Carlo Simulations
#' 
#' Calculates empirical p-value of t - or robust t - test
#' 
#' @author RT
#' @param x matrix of exogenous variables
#' @param y vector of endogenous variables
#' @param ymin lower bound of endogenous variable
#' @param ymax upper bound of endogenous variable
#' @param test specifies test either t - or robust t - test ("t" or "r" respectively)
#' @param alpha significance level of the test
#' @param betabrj supposedly true value of coefficient of interest under H_0
#' @param j index of the coefficient of interest
#' @param MC numeric value of Monte Carlo simulations
#' @param method defines what method is used to generate data: either division to the "Intervals" or k-nearest neighbors "KNN"
#' @param par in Intervals case defines in how many intervals range of Y is divided; in KNN case defines for how many nearest neighbours of each row to include in the generation process itself excluded.
#' @return empirical p-value
#' @examples
#' library(datasets)
#'  Z <- as.matrix(mtcars)
#'  X <- Z[,6]
#'  X <- cbind(rep(1, length(Z[,1])), X)
#'  Y <- Z[,1]
#'  #Calling method Intervals
#'  MCcont(X,Y, 10,40,"t",.05,0,2,10000,"Intervals",1) # for t-test, the range of Y is devided in 1 Interval
#'  MCcont(X,Y, 10,40,"r",.05,0,2,10000,"Intervals",3) # for robust t-test, the range of Y is devided in 3 Intervals
#'  ##Calling method KNN
#'  MCcont(X,Y, 10,40,"t",.05,0,2,10000,"KNN",3) # for t-test, with 3 nearest neighbors initialy assign to the generation process
#'  MCcont(X,Y, 10,40,"r",.05,0,2,10000,"KNN",6) # for robust t-test, with 6 nearest neighbors initialy assign to the generation process
#' @export
#' @useDynLib MCrt
MCcont<-function(x,y, ymin,ymax,test,alpha,betabarj,j,MC,method,par){
  #knn number of neighbors
  beta_hat<-solve(t(x)%*%x)%*%t(x)%*%y
  #rearnaging X and Y in a way that data are ordered accordingly to the values of a vector X*beta_hat increasingly
  ord<-order(X%*%beta_hat)
  X<-X[ord,]
  Y<-Y[ord]
  ## end of reordering
  tstat<-cppMCcont15(x,y,par,test,betabarj,j,MC, ymin,ymax,method)
  
  XROWS<-nrow(x)
  XCOLS<-ncol(x)
  if(length(tstat)<MC){
    cat(paste("Ups.. An error occured: There was at least one row for which it was unable to find more than one similar and feasible row with given parameters"))
  }
  else{
  #calculates emprical p-value from the t-statistics generated in C++ code
  emp_p<-sum(abs(tstat)>=qt(1-(alpha/2), df=(XROWS-XCOLS)))/MC 
    
  #cat(paste("Empirical p-value calculated with ",MC," MC simulations: ", "\n", sep=""))
  #cat(paste("                                                 t-test ",emp_p[1], "\n", sep=""))
  #cat(paste("                                                 robust t-test ",emp_p[2], "\n", sep=""))
  testw<-switch(test, t="t-test",r="robust t-test")
  Met<-switch(method, KNN="Data was generated using k-nearest neighbors method",Intervals="Data was generated using the division into the equal intervals of Y method")
  cat(paste("Empirical p-value for ",testw," calculated with ",MC," MC simulations : ", emp_p[1], "\n", sep=""))
  cat(paste(Met,"\n"))
  emp_p}
}

