#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[ Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

mat sigC3(mat dLtil, int K, int T, int N, vec id, vec H) {

  int i;
  int j;
  int t;
  int s;
  double temp1;
  double temp2;
  int Hi;
  int Hj;
  mat sig(K,K,fill::zeros);
  mat I(K,K,fill::eye);
  mat gi(K,1,fill::zeros);
  mat gj(K,1,fill::zeros);

  for (i=0;i!=N;i++)
    {
      for (j=0;j!=N;j++)
        {
          for (t=0;t!=T;t++)
            {
              for (s=0;s!=T;s++)
                {
                  temp1 = as_scalar(dLtil(i,t))*as_scalar(dLtil(i,s));
                  temp2 = (abs(t-s))/(T^(1/3));
  
                  if (fabs(temp2) <= 1)
                    {
                      temp2 = 1-fabs(temp2);
                    }
                  else
                    {
                      temp2 = 0;
                    }
                    
                  gi = I.col(id(i)-1);
                  gj = I.col(id(j)-1);
                  
                  temp1 = temp1*temp2;
                  Hi = H(id(i)-1);
                  Hj = H(id(j)-1);
                  temp2 = (Hi*Hj)^(1/2);
                  
                  sig = sig+((temp1/temp2)*(gi*gj.t())); 
                }
            }
        }
    }
    
  sig = sig/T;

  return sig;

}
