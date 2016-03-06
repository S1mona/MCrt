
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <string>
#include <vector>
#include <random>
using namespace Rcpp;
using namespace std;
using namespace arma;


// [[Rcpp::export]]
double funk1( arma::mat x, arma::vec y, String test,double betasbarj, int j)
{
    using namespace arma;
    using namespace std;
    int n = x.n_rows;
    int m = x.n_cols;
    
    arma::mat betas(m,1);
    arma::mat stand(1,1);
    //double stand;
    arma::vec t_stat(1);
    arma::mat resid(n,1);
    arma::mat xtx(m,m);
    
    xtx=(trans(x)*x).i();
    betas = xtx*trans(x)*y;
    resid=y-x*betas;
    if(test=="t")
    {
        stand=sqrt(xtx(j,j)*trans(resid)*resid/(n-m));
    }
    if(test=="r")
    {
        arma::mat omega(n,n);
        arma::mat stand1(m,m);
        omega=omega.fill(0);
        omega.diag()=square(resid);
        stand1=xtx*trans(x)*omega*x*xtx;
        stand=sqrt(stand1(j,j)*n/(n-m));
    }
    //small changes from original function smallest positive double is 2.2e-308 therefore we divide by it in else case.
    if(stand(0,0)>0)
    {
        t_stat=(betas(j)-betasbarj)/stand;
    }
    else
    {
        double dblmin = std::numeric_limits< double >::min();
        t_stat=(betas(j)-betasbarj)/dblmin;
    }
    return (double(t_stat(0)));
}


// [[Rcpp::export]]
NumericVector cppselect(arma::vec X,arma::vec index ) 
{ NumericVector resu(index.size());
  int i(0);
  for(i=0;i<index.size();i++)
  {resu(i)=X(index(i));
  }
return resu;
}


// [[Rcpp::export]]
NumericVector whichis(arma::vec X, double z) 
{
  int i(0);
  int b(0);
  arma::vec COcoef(X.n_elem);
  for(i=0;i<X.n_elem;i++){
    if(X(i)==z)
    {
       COcoef(i)=i;
       b++;
    }
   if(X(i)!=z)
   {
      COcoef(i)=-1;
    }
  }
    NumericVector resi(b);
    b=0;
    for(i=0;i<X.n_elem;i++)
    {
    if(COcoef(i)>0)
    {
       resi(b)=COcoef(i);
       b++;
    }
  
  }
  return resi;
}

// [[Rcpp::export]]
NumericVector inwhich_group(arma::mat Intervals,arma::vec v )
{
  int i(0);
  int b(0);
  NumericVector vv(v.n_elem);
  for(i=0;i<v.n_elem;i++){
      for(b=0;b<Intervals.n_rows;b++){
    if( (Intervals(b,0)<=v(i)) && (v(i)<=Intervals(b,1)) )
    {
       vv(i)=b;
       //i++;
    }
  
  }
  
  }
return vv;
}  

// [[Rcpp::export]]
double randomizer(NumericVector vvec) {
    //std::vector<int> numbers { 11 , 88 , -5 , 13 , 4 , 121 , 77 , 2 } ;
    std::random_device seed ;
    // generator
    std::mt19937 engine( seed( ) ) ;
    // number distribution
    std::uniform_int_distribution<int> choose( 0 , vvec.size( ) - 1 ) ;
 
 
    return vvec[ choose( engine ) ] ;
}




// [[Rcpp::export]]
arma::mat Sdefine(arma::vec xbetas, int intN,int n,arma::vec XBETASG, arma::vec resid, double ymin, double ymax){
    mat  S = zeros<mat>(n,n);
    int i(0);
    int b(0);
  //b) iv)
    // find out to which interval does xib belong
     // step 1 create the matrix with interval end values
    //double xbmin=xbetas.min();
    //double xbmax=xbetas.max();
    double intLength=(ymax-ymin)/intN;
    arma::mat Inter(intN,2);
    // loop creates matrix with the end points of the intervals. It is used to groups rows to similar rows
    for(i=0;i<intN;i++)
    {
        Inter(i,0) = ymin+intLength*i;
        Inter(i,1)= ymin+intLength*(i+1);
    }
 
    // Information for which group each row (aka X*betas) belongs is stored in vector gr
    NumericVector gr=inwhich_group(Inter,xbetas);
   // fill S matrix with ones in (i,j) if i and j are in the same group
   
   for(i=0;i<n;i++)
    { 
      NumericVector indi=whichis(gr,gr(i));
      for(b=0;b<indi.size();b++)
      {
        S(i,indi(b))=1;
      }
    }
   
    /* Here we run a loop over all rows and check if all residuals of similar rows are feasible (between 
   the bound of Y) for that row. If one or few residudals are not feasible, not feasible residual with the
   largest absolute value is taken from the group. And new group with new mean is tested for the same,
   untill all residual for row i are similar and feasible.
   */

    double yy; //initiation of variable Y used in while loop to calculate feasibility of residuals
    int calc(0); //initiation of variable calc which is used in while loop to check if all smilar rows
                    //are in bounds of Y
    double M_rfori; // mean of feasible and in group residuals for i-th row
    
    for(i=0;i<n;i++)
    {// the loop runs through all rows
        calc=0;
        NumericVector indii=whichis(trans(S.row(i)),1); //vector of indices of similar rows
        NumericVector rfori=cppselect(resid,indii); //select actual residuals from similar group
        M_rfori=mean(rfori);// mean of similar residuals
       
       while (calc<indii.size()) {//while loop runs until all similar rows are feasible
            int prob_res_ind=-1; //initiation of variable. Here index of row which is 
                                 //not feasible and has absolute largest residual is stored
            calc=0;
            
            for(b=0;b<indii.size();b++)
            {//run through all similar residuals and if it is feasible add 1 to calc
                
              yy=XBETASG(i)+resid(indii(b))-M_rfori;
                
                if (yy<=ymax && yy>=ymin ) {
                    calc=calc+1;
                }
                else
                { 
                      if(prob_res_ind==-1)
                      {
                          prob_res_ind=indii(b);
                      }
                      else
                      {
                        if(abs(resid(prob_res_ind))<abs(resid(indii(b)))){
                          prob_res_ind=indii(b);
                        }
                      }
                }
            }
           
           if ( prob_res_ind > -1) {
                S(i,prob_res_ind)=0;
                indii=whichis(trans(S.row(i)),1);
                rfori=cppselect(resid,indii); //select actual residuals from similar group
                M_rfori=mean(rfori);// mean of similar residuals 
            }
      } //end of while loop
    }// Feasibility is checked for all rows.
   return S;
}

// defining S matrix for knn case
// [[Rcpp::export]]
arma::mat Sdefine_knn(arma::vec xbetas, int knn,int n,arma::vec XBETASG, arma::vec resid, double ymin, double ymax){
    mat  S = eye<mat>(n,n);
    int i(0);
    int b(0);
    int d(0);
    
    double dif1;
    double dif2;
    arma::vec nei(2);
   // fill S matrix with ones in (i,j) if i and j are in the same group
   
   for(i=0;i<n;i++){
     std::fill(nei.begin(), nei.end(), 0);// setting nei vector to zeros for new row
     nei(0)=0;
     nei(1)=0;
     for(b=1;b<=knn;b++){
       if(i-(nei(0)+1)<0){nei(1)=nei(1)+1;}// checking if we do not exied the boundaries of indices from below
        else{
            if(i+(nei(1)+1)>=n){nei(0)=nei(0)+1;}// checking if we do not exied the boundaries of indices from above
              else{//looking for the nearest neighbour in regard to xbetas value with priority to rhs
                  dif1=abs(xbetas(i)-xbetas(i-(nei(0)+1)));
                  dif2=abs(xbetas(i+(nei(1)+1))-xbetas(i));
                  if(dif1>dif2){nei(1)=nei(1)+1;}
                  else{nei(0)=nei(0)+1;}
       
                  }
            }
     }
    for(d=0;d<=nei(0);d++){ 
       S(i,i-d)=1;
     }//assigning ones in the indices Matrix if other rows are in the neigberhoud on the lhs
     
     for(d=0;d<=nei(1);d++){
       S(i,i+d)=1;
     } //assigning ones in the indices Matrix if other rows are in the neigberhoud on the rhs
   
   }
   
    /* Here we run a loop over all rows and check if all residuals of similar rows are feasible (between 
   the bound of Y) for that row. If one or few residudals are not feasible, not feasible residual with the
   largest absolute value is taken from the group. And new group with new mean is tested for the same,
   untill all residual for row i are similar and feasible.
   */

    double yy; //initiation of variable Y used in while loop to calculate feasibility of residuals
    int calc(0); //initiation of variable calc which is used in while loop to check if all smilar rows
                    //are in bounds of Y
    double M_rfori; // mean of feasible and in group residuals for i-th row
    
    for(i=0;i<n;i++)
    {// the loop runs through all rows
        calc=0;
        NumericVector indii=whichis(trans(S.row(i)),1); //vector of indices of similar rows
        NumericVector rfori=cppselect(resid,indii); //select actual residuals from similar group
        M_rfori=mean(rfori);// mean of similar residuals
       
       while (calc<indii.size()) {//while loop runs until all similar rows are feasible
            int prob_res_ind=-1; //initiation of variable. Here index of row which is 
                                 //not feasible and has absolute largest residual is stored
            calc=0;
            
            for(b=0;b<indii.size();b++)
            {//run through all similar residuals and if it is feasible add 1 to calc
                
              yy=XBETASG(i)+resid(indii(b))-M_rfori;
                
                if (yy<=ymax && yy>=ymin ) {
                    calc=calc+1;
                }
                else
                { 
                      if(prob_res_ind==-1)
                      {
                          prob_res_ind=indii(b);
                      }
                      else
                      {
                        if(abs(resid(prob_res_ind))<abs(resid(indii(b)))){
                          prob_res_ind=indii(b);
                        }
                      }
                }
            }
           
           if ( prob_res_ind > -1) {
                S(i,prob_res_ind)=0;
                indii=whichis(trans(S.row(i)),1);
                rfori=cppselect(resid,indii); //select actual residuals from similar group
                M_rfori=mean(rfori);// mean of similar residuals 
            }
      } //end of while loop
    }// Feasibility is checked for all rows.
   return S;
}

//' @export
// [[Rcpp::export]]
arma::mat cppMCcont15(arma::mat x, arma::vec y, int intN, String test, double betabarj, int j, int MC, double ymin, double ymax,String method){
    int knn = intN;
    //initiating global function variables
    int n = x.n_rows; // number of  data rows / size N
    arma::mat t_stat_all(1,MC); // vector in which t-stats are saved
    mat  S = zeros<mat>(n,n);
    int i(0);
    int b(0);
    
    // a) estimating beta hat & residuals ols/robust case
    arma::mat xtx=(trans(x)*x).i();
    arma::mat betas = xtx*trans(x)*y;
    arma::vec xbetas=x*betas;
    arma::vec resid=y-xbetas;
    
    arma::vec betasg=betas; // beta vector used to simulate Ys
    betasg.at(j-1)=betabarj; // j-th covariate is assigned to the true value
    arma::vec XBETASG=x*betasg; //deterministic value of x*betasg also used to generate Ys
    arma::vec ErMean(n); //initiation of ErMean for each group
    
    
    // S matrix is defined in the line below. The ones in particular row of this matrix indicate other similar rows with feasible residuals to it.
    if(method=="Intervals"){
    S=Sdefine(xbetas, intN,n, XBETASG, resid, ymin, ymax);}
    if(method=="KNN"){
    S=Sdefine_knn(xbetas, knn,n,XBETASG, resid,  ymin,  ymax);}
    
    //run through all S rows and calculate ErMean
    for(i=0;i<n;i++){
      ErMean(i)=mean(cppselect(resid,whichis(trans(S.row(i)),1)));// select actual residuals from similar group and calc mean of similar residuals
       // mean value of residuals from similar and feasible rows  for i-th row
    }
    //check if all rows have at least one similar row(not equal to itself)
    arma::vec checkR(n);
    for (i=0; i<n; i++) {
        checkR(i)=sum(S.row(i));
    }
    if(checkR.min()>=2){
        ////////////////
        // generate new Y population with beta_J equal to the true value
        ///////////////
        arma::vec RandErr(n); // initiate variable for random errors
        arma::vec new_Y(n); // initiate variable for  new Y realisation
        
        for(b=0;b<MC;b++){
            //MC simulations loop
            
            for(i=0;i<n;i++)
            {//generate vector of random Error for each
                NumericVector simRow=whichis(trans(S.row(i)),1);
                RandErr(i)= resid(randomizer(simRow));
            }
            
            
            // generate new Y realisation:
            new_Y=x*betasg+RandErr-ErMean;
            //calculate t statistic for j-th covariate  call function funk1
            t_stat_all(0,b)= funk1(x,new_Y,test,betabarj,(j-1)) ;
            //t_stat_all(1,b)= funk1(x,new_Y,"r",betabarj,(j-1)) ;
            
        }//end of MC loop
        return t_stat_all;
    }
    else{
        arma::mat error_code(1,1);
        error_code(0,0)=101;
        return error_code; //error code for too few similar rows, ask user to decrease the number of intervals
    }
  //return S;
}