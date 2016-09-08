#include <stdio.h>
#include <math.h>
#include <complex.h>
//testlib.c

double choose(int n, int k){
  if (k>n && k<0){
    return 0;
  }
  else if (k > n/2){
    return choose(n,n-k);
  }
  else if (k == 0){
    return 1;
  }
  else{
    return n*choose(n-1,k-1)/k;
  }    
}

double complex e(int N, int k, double x, double y){
  int i;
  double complex result = sqrt(choose(N,k));

  result /= pow(1.0+x*x+y*y,(double)N/2.0);
  for(i=0; i<k; i++){
    result *= x+y*I;
  }
  return result;
}

double h(double x, double y){
  return pow(x*x+y*y-1.0,2.0)*(1.0+x*x+y*y+x)/pow(x*x+y*y+1.0,3.0);
  //return 1.0;
}

double realIntegrand(int N, int k, int l, double x, double y){
  return creal(e(N,k,x,y)*h(x,y)*conj(e(N,l,x,y))/pow(1.0+x*x+y*y,2.0));
}

double imagIntegrand(int N, int k, int l, double x, double y){
  return cimag(e(N,k,x,y)*h(x,y)*conj(e(N,l,x,y))/pow(1.0+x*x+y*y,2.0));
}

