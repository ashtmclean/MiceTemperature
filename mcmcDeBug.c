#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include <stdbool.h>

#define VSET gsl_vector_set
#define VGET gsl_vector_get
#define MSET gsl_matrix_set
#define MGET gsl_matrix_get

double f1t(double theta, double B, int t){
  double val = (1/(1 + pow((12/M_PI),2)*pow(theta,2)))*(12/M_PI)*(sin(((M_PI*t)/12)+B)+(12/M_PI)*theta*cos(((M_PI*t)/12)+B));
  return(val);
}


double f2t(double theta, double B, int t){
  double val = (1/(1 + pow((12/M_PI),2)*pow(theta,2)))*(12/M_PI)*(sin(((M_PI*t)/12)+B)-(12/M_PI)*theta*cos(((M_PI*t)/12)+B));
  return(val);
}

double g0t(double theta, double B, int t){
  double s = sin((M_PI*t)/12 + B);
  double s2 = sin(2*((M_PI*t)/12 + B));
  double v = pow((12/M_PI),2)*(1/(1 +pow((12/M_PI),2)*pow(theta,2)))*(0.5*pow(s,2) + theta*((t + (12/(2*M_PI))*s2)/2));
  return(v);
}


double g1t(double theta, double B, int t){
  double tpi = 12/(2*M_PI);
  double targ = 2*(((M_PI*t)/12) + B);
  double val = (1/(1 + (4*(pow(tpi,2)))*(pow(theta,2))))*tpi*(-0.5*cos(targ) - tpi*theta*sin(targ));
  return(val);
}

double g2t(double theta, double B, int t){
  double tpi = 12/(2*M_PI);
  double targ = 2*(((M_PI*t)/12) + B);
  double val = (1/(1 + (4*pow(tpi,2))*pow(theta,2)))*tpi*(0.5*sin(targ) - tpi*theta*cos(targ));
  return(val);
}

double g3t(double theta, double B, int t){
  double tpi = 12/M_PI;
  double coef = 1/(1 + pow(tpi,2)*pow(theta,2));
  double ex = exp(-2*theta*t);
  double val = coef*tpi*(g1t(theta, B, t) - tpi*theta*((-1/(4*theta))*ex + g2t(theta, B, t)));
  return(val);
}



double v11t(gsl_vector *param, gsl_matrix *v0, int t, int tm){
  double theta1 = gsl_vector_get(param, 0);
  double sig = gsl_vector_get(param, 2);
  double sigA = gsl_vector_get(param, 3);
  double B = gsl_vector_get(param, 4);
  
  double v1 = gsl_matrix_get(v0, 0, 0);
  double v2 = gsl_matrix_get(v0, 0, 1);
  double v3 = gsl_matrix_get(v0, 1, 1);
  
  double ic = 2*v2*(M_PI/12)*(exp(-theta1*(t-tm))*f2t(theta1, B, t)-f2t(theta1, B, tm));
  double coef = pow((M_PI/12),2)*(pow(sigA,2)/theta1);
  double p4 = -1/(4*theta1)*(exp(-2*theta1*(t-tm))-1) + exp(-2*theta1*(t-tm))*g2t(theta1, B, t) - g2t(theta1, B, tm);
double p1 = g0t(theta1, B, t) - g0t(theta1, B, tm) - exp(-theta1*(t-tm))*f1t(theta1, B, tm)*f2t(theta1, B, t)+f1t(theta1, B, tm)*f2t(theta1, B, tm);

 double p2 = exp(-theta1*(t-tm))*f2t(theta1, B, tm)*f2t(theta1, B, t)-f2t(theta1, B, tm)*f2t(theta1, B, tm);

 double p3 = (1/(1+ (pow((12/M_PI),2)*pow(theta1,2))))*(12/M_PI)*(exp(-2*theta1*(t-tm))*g1t(theta1, B, t) - g1t(theta1, B, tm)-(12/M_PI)*theta1*p4);
  double ic2 = 2*v3*pow((M_PI/12),2)*(p3-p2);
  double v = v1+ ic + ic2 + coef*(p1 + p2 - p3) + pow(sig,2)*(t - tm);
 //printf("\n%f\n", v);
  return(v);
}

double v22t(gsl_vector *param, gsl_matrix *v0, int t, int tm) 
{
  double theta1 = gsl_vector_get(param, 0);
  double sigA = gsl_vector_get(param, 3);  
  double v3 = gsl_matrix_get(v0, 1, 1);
  double v = v3*exp(-2*theta1*(t-tm))+(pow(sigA,2)/(2*theta1))*(1-exp(-2*theta1*(t-tm)));  
  
  return (v);
}

double v12t(gsl_vector *param, gsl_matrix *v0, int t, int tm)
{
  double theta1 = gsl_vector_get(param, 0);
  double sigA = gsl_vector_get(param, 3);
  double B = gsl_vector_get(param, 4); 
  double v1 = gsl_matrix_get(v0, 0, 0);
  double v2 = gsl_matrix_get(v0, 0, 1);
  double v3 = gsl_matrix_get(v0, 1, 1);
  double new=(M_PI/12)*v3*exp(theta1*(2*(tm-t)))*f2t(theta1,B,t)-(M_PI/12)*v3*exp(theta1*(tm-t))*f2t(theta1,B,tm);
  double v = v2*exp(-theta1*(t-tm))+new +(f1t(theta1,B,t)-exp(-theta1*(t-tm))*f1t(theta1,B,tm)-exp(-theta1*2*(t-tm))*f2t(theta1,B,t)+exp(-theta1*(t-tm))*f2t(theta1,B,tm))*((M_PI*pow(sigA,2))/(24*theta1));
  return (v); 
}

double m2t(gsl_vector *param, gsl_vector *m0, int t, int tm){
  double theta1 = gsl_vector_get(param, 0);
  double theta2 = gsl_vector_get(param, 1); 
  double init=gsl_vector_get(m0, 1);
  double v = init*exp(-theta1*(t-tm))+theta2 * (1 - exp(-theta1*(t-tm)));
  return (v); 
}

double m1t(gsl_vector *param, gsl_vector *m0, int t,int tm){
  double theta1 = gsl_vector_get(param, 0);
  double theta2 = gsl_vector_get(param, 1);
  double B = gsl_vector_get(param, 4);	
  double init1 =gsl_vector_get(m0, 0);
  double init2 =gsl_vector_get(m0, 1);
  double v = init1+(M_PI/12)*(theta2*(12/M_PI)*sin(((M_PI*t)/12)+B) - theta2*(12/M_PI)*sin(((M_PI*tm)/12)+B) + (init2-theta2)*(f2t(theta1,B,t)*exp(theta1*(tm  - t)) -f2t(theta1,B,tm)));//exp(theta1*tm)*(init2-theta2)*(f2t(theta1,B,t) -f2t(theta1,B,tm)));
  
 //(m0[2]-theta2)*(f2tb(theta1,B,t)*exp(theta1*(tm  - t)) -f2tb(theta1,B,tm)))
  
  return(v);
}

double lnd(double x,double mean,double sd)
{
  /* returns the log of a gaussian density with "mean" and "sd" */
  
  double a;
  a=-0.5*log(2*M_PI)-log(sd)-0.5*((x-mean)*(x-mean)/(sd*sd));
  return(a);
}

double lmgpdf(int d,gsl_vector *x,gsl_vector *mu,gsl_matrix *var)
{
  /* Returns the log of a multivariate gaussian pdf with mean vec, mu and var matrix, var */
  /* up to an additive constant */

  int i;
  double det=1.0;
  double s=0.0;
  double ll;
  gsl_matrix *disp; 
  gsl_vector *temp;
 
  disp=gsl_matrix_alloc(d,d);
  temp=gsl_vector_alloc(d);

  gsl_matrix_memcpy(disp,var);
  gsl_linalg_cholesky_decomp(disp);
  
  for (i=0;i<d;i++)
    {
      det=det*MGET(disp,i,i);
      VSET(temp,i,VGET(x,i)-VGET(mu,i));
    }

  gsl_linalg_cholesky_svx(disp,temp);

  for (i=0;i<d;i++)
    {
      s=s+(VGET(temp,i))*(VGET(x,i)-VGET(mu,i));
    }
  ll=-1.0*log(det)-0.5*s;

  gsl_vector_free(temp);
  gsl_matrix_free(disp);

  return(ll);
}



double logPrior(gsl_vector *param){
	double lamb1 = (gsl_vector_get(param, 0));
	double lamb2 = (gsl_vector_get(param, 1));
	double lamb3 = exp(gsl_vector_get(param, 2));
	double lamb4 = exp(gsl_vector_get(param, 3));
	double lamb5 = (gsl_vector_get(param, 4));
	double lamb6 = exp(gsl_vector_get(param, 5));
	double sumLamb = lamb3 + lamb4 +lamb6;	
	double ll = (lnd(lamb1, -1, 1)) + (lnd(lamb2, 1.2, 0.3)) + log(pow(lamb3, -4)*exp(-(1/lamb3))) + log(pow(lamb4, -4)*exp(-(1/lamb4)))+ log(gsl_ran_flat_pdf(lamb5, -(M_PI), (M_PI)))+ log(pow(lamb6, -4)*exp(-(3/lamb6))) + sumLamb;;  
 	return ll;
}

double dotProd(gsl_vector *a, gsl_vector *b){
	int i;
	int al = a->size;
	int bl = b->size;
	double tot = 0;
	if(al == bl){
		for (i=0; i< al; i++){
			tot = tot + gsl_vector_get(a, i)*gsl_vector_get(b, i);
		}
		return tot;
	}
	return 0;
}

void matmul(int d1,int d2,int d3,gsl_matrix *A,gsl_matrix *B,gsl_matrix *X)
{
  /* multiplies A*B where A is d1*d2 and B is d2*d3 - places result in C */

  int i,j,k;
  
  //gsl_matrix_set_all(C,0.0);
  for (i=0;i<d1;i++)
    {
      for(j=0;j<d3;j++)
		{
 			for(k=0;k<d2;k++)
  			 {
     			MSET(X,i,j,MGET(X,i,j)+MGET(A,i,k)*MGET(B,k,j));
   				}
			}
    }
}



void vecmatmul(int d1, int d2, gsl_vector *B, gsl_matrix *A, gsl_vector *C)
{
  /* multiplies A*B where A is d1*d2 and B is d2*1 - places result in vector C */

  	int i,k;
  	//printf("\n%d\n", d1);
 	gsl_vector_set(C, 0, 0.0);
 	gsl_vector_set(C, 1, 0.0);

  	for (i=0;i<d1;i++)
    {
      for(k=0;k<d2;k++)
		{
		gsl_vector_set(C, i,gsl_vector_get(C, i) + gsl_matrix_get(A, k, i)*gsl_vector_get(B, k));
 		//VSET(C,i,VGET(C,i)+MGET(A,i,k)*VGET(B,k));
		}
    }
}

void invert(int d,gsl_matrix *A,gsl_matrix *B)
{ 
  int i;
  double det; 
  gsl_vector *s;
  gsl_vector *work;
  gsl_matrix *temporary;
  gsl_matrix *temporary2; 
  gsl_matrix *V;
  s=gsl_vector_alloc(d);
  work=gsl_vector_alloc(d);
  temporary=gsl_matrix_alloc(d,d);
  temporary2=gsl_matrix_alloc(d,d);
  V=gsl_matrix_alloc(d,d); 
  
  /* Inverts A and places result in B */ 
  if(d==2){
    det = MGET(A,0,0)*MGET(A,1,1)-MGET(A,0,1)*MGET(A,1,0);
    MSET(B,0,0,MGET(A,1,1)/det);
    MSET(B,1,1,MGET(A,0,0)/det);
    MSET(B,0,1,-1.0*MGET(A,0,1)/det);
    MSET(B,1,0,-1.0*MGET(A,1,0)/det);
  }else if(d==1){
    MSET(B,0,0,1.0/(MGET(A,0,0)));
  }else{
  gsl_matrix_memcpy(B,A);
  gsl_linalg_SV_decomp(B,V,s,work);
  
  gsl_matrix_set_all(temporary,0.0);
  gsl_matrix_set_all(temporary2,0.0);
  for(i=0;i<d;i++)
    {
      MSET(temporary,i,i,(1.0/VGET(s,i)));
    }
  gsl_matrix_transpose(B);
  matmul(d,d,d,temporary,B,temporary2);
  matmul(d,d,d,V,temporary2,B);
  }
  
  gsl_vector_free(s);
  gsl_vector_free(work);
  gsl_matrix_free(temporary);
  gsl_matrix_free(temporary2); 
  gsl_matrix_free(V);
}

void matvecmul(int d1, int d2,  gsl_matrix *A, gsl_vector *B, gsl_vector *C)
{
  /* multiplies A*B where A is d1*d2 and B is d2*1 - places result in vector C */

  	int i,k;
  	//printf("\n%d\n", d1);
 	gsl_vector_set(C, 0, 0.0);
 	gsl_vector_set(C, 1, 0.0);
  	for (i=0;i<d1;i++)
    {
      for(k=0;k<d2;k++)
		{
		gsl_vector_set(C, i,gsl_vector_get(C, i) + gsl_matrix_get(A, i, k)*gsl_vector_get(B, k));
 		//VSET(C,i,VGET(C,i)+MGET(A,i,k)*VGET(B,k));
		}
    }
}

void vecScale(gsl_vector *a, double b, gsl_vector *ab){
	int i; 
	int j = a -> size;
	for(i=0; i<j; i++){
	gsl_vector_set(ab, i, gsl_vector_get(a, i)*b);
	}
}

void vecTvec(gsl_vector *a, gsl_vector *b, gsl_matrix *c){
	int i, j;
	int n = a -> size; 
	int m = b -> size;
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			gsl_matrix_set(c, i, j, gsl_vector_get(a, i)*gsl_vector_get(b, j));
		}
	}

}

void mymatmul(gsl_matrix *C1, gsl_matrix *Vt, gsl_matrix *C2){
	double a, b, c, d;
	gsl_vector *Cab, *Ccd, *Vac, *Vbd;
	Cab = gsl_vector_alloc(2);
	Ccd = gsl_vector_alloc(2);
	Vac = gsl_vector_alloc(2);
	Vbd = gsl_vector_alloc(2);
	
	gsl_vector_set(Cab, 0, gsl_matrix_get(C1, 0, 0));
	gsl_vector_set(Cab, 1, gsl_matrix_get(C1, 0, 1));
	
	gsl_vector_set(Ccd, 0, gsl_matrix_get(C1, 1, 0));
	gsl_vector_set(Ccd, 1, gsl_matrix_get(C1, 1, 1));
	
	gsl_vector_set(Vac, 0, gsl_matrix_get(C1, 0, 0));
	gsl_vector_set(Vac, 1, gsl_matrix_get(C1, 1, 0));
	
	gsl_vector_set(Vbd, 0, gsl_matrix_get(C1, 0, 1));
	gsl_vector_set(Vbd, 1, gsl_matrix_get(C1, 1, 1));
	
    a = dotProd(Cab, Vac);
  	b = dotProd(Cab, Vbd);
  	c = dotProd(Ccd, Vac);
  	d = dotProd(Ccd, Vbd);
  	
	gsl_matrix_set(C2, 0,0, a);
 	gsl_matrix_set(C2, 0,1, b);
 	gsl_matrix_set(C2, 1,0, c);
 	gsl_matrix_set(C2, 1,1, d);
 	
 	gsl_vector_free(Cab);
 	gsl_vector_free(Ccd);
 	gsl_vector_free(Vac);
 	gsl_vector_free(Vbd);
}

double filter(gsl_vector *param, gsl_vector *a, gsl_matrix *C, gsl_vector *F, gsl_matrix *X){

	
	gsl_vector_set(a,0, 36);
    gsl_vector_set(a,1, 5);

    gsl_matrix_set(C, 0, 0, 1);
    gsl_matrix_set(C, 0, 1, 0);
    gsl_matrix_set(C, 1, 0, 0);
    gsl_matrix_set(C, 1, 1, 0.1);
    
    gsl_vector_set(F, 0, 1);
	gsl_vector_set(F, 1, 0);

	double lamb1 = gsl_vector_get(param, 0);
	double lamb2 = gsl_vector_get(param, 1);
	double lamb3 = gsl_vector_get(param, 2);
	double lamb4 = gsl_vector_get(param, 3);
	double lamb5 = gsl_vector_get(param, 4);
	double lamb6 = gsl_vector_get(param, 5);
	int n = X->size1; // number of rows (observations)

	gsl_vector *mt; 
	mt = gsl_vector_alloc(2);
	gsl_matrix *Vt;
	Vt = gsl_matrix_alloc(2,2);
	double ll = 0;
	double mean1 = dotProd(F, a);
	gsl_vector *FC;
	FC = gsl_vector_alloc(2);
	gsl_vector *CF;
	CF = gsl_vector_alloc(2);
	int d1 = C->size1;
	int d2 = C->size2;
	vecmatmul(d1, d2, F, C, FC);

	matvecmul(d1, d2, C, F, CF);
	double FCF = dotProd(FC, F);
	double sd = sqrt(FCF + lamb6);

	ll = ll+ lnd(gsl_matrix_get(X, 0,0), mean1, sd);

	//UPDATE POSTERIOR A
	gsl_vector *a1; 
	a1 = gsl_vector_alloc(2);
	vecScale(CF, (1/(FCF+pow(lamb6, 2))), a1);
	double a2;
	a2 = gsl_matrix_get(X, 0, 0) - dotProd(F, a);
	gsl_vector *a3; 
	a3 = gsl_vector_alloc(2);
	gsl_vector *a4; 
	a4 = gsl_vector_alloc(2);
	vecScale(a1, a2, a3);
	gsl_vector_add(a, a3); 

	// B
	gsl_matrix *C1;
	C1 = gsl_matrix_alloc(2,2);
	vecTvec(a1, F, C1);
	gsl_matrix *C2;
	C2 = gsl_matrix_alloc(2,2);	
	matmul(2, 2, 2, C1, C, C2);
	gsl_matrix_sub(C, C2);

	gsl_matrix *V2;
	V2 = gsl_matrix_alloc(2,2);	

	int i;
	for(i=0; i<(n-1); i++){
		int tm = gsl_matrix_get(X, i, 1);
		int t = gsl_matrix_get(X, (i+1), 1);

		//UPDATE m
		
		gsl_vector_set(mt, 0, m1t(param, a, t, tm));
		gsl_vector_set(mt, 1, m2t(param, a, t, tm));

		//UPDATE V
		gsl_matrix_set(Vt, 0, 0, v11t(param, C, t, tm));
		gsl_matrix_set(Vt, 1, 1, v22t(param, C, t, tm));
		gsl_matrix_set(Vt, 0, 1, v12t(param, C, t, tm));
		gsl_matrix_set(Vt, 1, 0, gsl_matrix_get(Vt, 0, 1));
		
		//printf("\n%d\n", i);
		//printf("\n%f\n", gsl_matrix_get(Vt, 1, 1));
		
		// UPDATE LL
		mean1 = dotProd(F, mt);

		vecmatmul(d1, d2, F, Vt, FC);
		FCF = dotProd(FC, F);
		sd = sqrt(FCF + lamb6);

		ll = ll + lnd(gsl_matrix_get(X, (i+1),0), mean1, sd);
		
		matvecmul(d1, d2, Vt, F, CF);
		vecmatmul(d1, d2, F, Vt, FC);
		double FCF = dotProd(FC, F);
		vecScale(CF, (1/(FCF+pow(lamb6, 2))), a1);
		a2 = gsl_matrix_get(X, (i+1), 0) - dotProd(F, mt);
		vecScale(a1, a2, a3);
		gsl_vector_memcpy(a4,mt);
		gsl_vector_add(a4, a3);
		gsl_vector_memcpy(a, a4);

		vecTvec(a1, F, C1);
		gsl_matrix_set_all(C2, 0);
		matmul(2, 2, 2, C1, Vt, C2);

		gsl_matrix_memcpy(V2, Vt);
		gsl_matrix_add(V2, C2);
		gsl_matrix_memcpy(C, V2);
		
		
		
	}


	gsl_vector_free(mt);
 	gsl_matrix_free(Vt);
	gsl_vector_free(FC);
	gsl_vector_free(CF);
  	gsl_vector_free(a1);
  	gsl_vector_free(a3);
	gsl_matrix_free(C1);
	gsl_matrix_free(C2);
	gsl_matrix_free(V2);
	//printf("\n%f\n", ll);
	return ll;
}

void rmvn(gsl_vector *vec, gsl_vector *candidate, gsl_matrix *var, int size, gsl_rng *r){
 		
  		int i, j;
  		gsl_matrix *disp;
  		gsl_vector *ran;
  		gsl_vector *x;
  		ran=gsl_vector_alloc(size);
  		disp = gsl_matrix_alloc(size, size);
  		gsl_matrix_memcpy(disp,var);
  		x=gsl_vector_alloc(size);
  		gsl_vector_set_all(x,0.0); // do i need this line
  		gsl_linalg_cholesky_decomp(disp);
  		// make it triangular
  		for (i=0;i<size;i++){
    		for (j=i+1;j<size;j++){
 				MSET(disp,i,j,0.0);
			}
    	}
    	
    	
  		for (j=0;j<size;j++){
    		gsl_vector_set(ran,j,gsl_ran_gaussian(r,1.0));
    	}
    	
    	for (i=0;i<size;i++){
     		for (j=0;j<size;j++){
 				gsl_vector_set(x,i,gsl_vector_get(x,i)+gsl_matrix_get(disp,i,j)*gsl_vector_get(ran,j));
			}
    	}
    	for(i=0;i<size;i++){
    	gsl_vector_set(candidate,i, gsl_vector_get(vec, i)+gsl_vector_get(x,i)); //add mean
  	}

  	
}

void mcmcGrid(int iters, gsl_matrix *X, gsl_vector *init, gsl_matrix *sigtune, gsl_vector *x0, gsl_matrix *C, gsl_vector *F){
	
	int rej = 0;
	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);   
	int n = X->size1; // number of rows (observations)
	int p = init->size; //number of parameters
	int i, j, k;
	double margllcurr, margllcan, laprob;
	gsl_matrix *mat;
	mat = gsl_matrix_alloc(iters, (p+1));
	gsl_vector *current;
	current = gsl_vector_alloc(p);
	gsl_vector *can;
	can = gsl_vector_alloc(p);
	gsl_vector *expcan;
	expcan = gsl_vector_alloc(p);
	for(i = 0; i<p; i++){
		if(i == 4){
			gsl_matrix_set(mat, 0, i, gsl_vector_get(init, i));
			gsl_vector_set (current, i, gsl_vector_get(init, i));
			
		}
		else{
			gsl_matrix_set(mat, 0, i, gsl_vector_get(init, i));
			gsl_vector_set (current, i, log(gsl_vector_get(init, i)));
			
		}
	
	
	gsl_matrix_set(mat, 0, 6, 0);
	margllcurr = -100000;


	for(i=1;i<iters; i++){
		rmvn(current, can, sigtune, p, r);
		
	
		if(gsl_vector_get(can, 4)>(M_PI)){
			gsl_vector_set(can, 4, (-2*M_PI + gsl_vector_get(can, 4)));
			
		}
		
		if(gsl_vector_get(can, 4)<(-M_PI)){
			gsl_vector_set(can, 4, (2*M_PI + gsl_vector_get(can, 4)));
		}
		
		for(j=0;j<p; j++){
			if(j == 4){
				gsl_vector_set(expcan, j, (gsl_vector_get( can, j)));
			}
			else{
				gsl_vector_set(expcan, j, exp(gsl_vector_get( can, j)));
			}
		}
		
		margllcan = filter(expcan, x0, C, F, X);
			
	// DEBUGGING WHY HARD CUT OFF FOR THETA 1 
	
	//double a;
	//a = exp(gsl_vector_get(can, 0));
	
	//if(a > 0.6){
		//printf("\n%f %f %f %f %f %f %f", gsl_vector_get(expcan, 0), gsl_vector_get(expcan, 1), gsl_vector_get(expcan, 2), gsl_vector_get(expcan, 3), gsl_vector_get(expcan, 4), gsl_vector_get(expcan, 5), margllcan);
	//}
	
	// END DEBUGGING
		if(margllcan != margllcan){
			margllcan = -pow(10, 10);
		}
		//printf("\n%f\n", margllcan);

	laprob = margllcan - margllcurr + logPrior(can) - logPrior(current);
	
	

	
	if(log(gsl_ran_flat(r, 0, 1))< laprob){
		margllcurr = margllcan;
		gsl_vector_memcpy (current, can);
	}
	
	for(k=0;k<(p+1);k++){
		if(k==4){
			gsl_matrix_set(mat, i, k, gsl_vector_get(current, k));
		}
		else if(k==6){
			gsl_matrix_set(mat, i, k, margllcurr);
		}
		else{
			gsl_matrix_set(mat, i, k, exp(gsl_vector_get(current, k)));
		}
	}
	
	 //after debugging, uncomment the line below
	
		 printf("\n%f %f %f %f %f %f %f", gsl_matrix_get(mat, i,0), gsl_matrix_get(mat, i,1), gsl_matrix_get(mat, i,2), gsl_matrix_get(mat, i,3), gsl_matrix_get(mat, i,4), gsl_matrix_get(mat, i,5), gsl_matrix_get(mat, i,6));
	
		
	}
	
	

	}


}


int main( int argc, char *argv[] ){
	
	char *File1, *File2;
	int arg1 = atoi(argv[1]);
	int arg3 = atoi(argv[3]);
	int iters, i;
	//printf("/n%s/n", "hello");
	gsl_vector *init, *a, *F;
	gsl_matrix *data, *C;
	gsl_matrix *sigtune;
	sigtune = gsl_matrix_alloc(6,6);
	
	if(arg1 == 1){
	
		iters = 10000;
		gsl_matrix_set_zero (sigtune);
		for(i=0; i<6; i++){
			gsl_matrix_set(sigtune, i, i, 0.01);
		
		}
	}
	
	if(arg1 == 2){
		
		iters = 100000;
		FILE *pv;
		File2="postvar.dat"; 
  		pv=fopen(File2,"r"); 	
  		
 		gsl_matrix_fscanf(pv,sigtune);
 		fclose(pv);
 		
	}
	
	init=gsl_vector_alloc(6);

	gsl_vector_set(init,0, 0.1);
    gsl_vector_set(init,1, 3);
    gsl_vector_set(init,2, 0.1);
    gsl_vector_set(init,3, 0.1);
    gsl_vector_set(init,4, 0);    
    gsl_vector_set(init,5, 1);

    data = gsl_matrix_alloc(arg3, 2);

	F=gsl_vector_alloc(2);
	a=gsl_vector_alloc(2);
	C=gsl_matrix_alloc(2,2);
	
  	FILE *f; 
  	//Read in the data
  	f = fopen(argv[2], "r");
 	gsl_matrix_fscanf(f,data);
 	
 	fclose(f); 
	
	//filter(init, a, C, F, data);
	mcmcGrid(iters, data, init, sigtune, a, C, F);
 	gsl_matrix_free(data);
	gsl_matrix_free(C);
	gsl_vector_free(init);
  	gsl_vector_free(a);
 	gsl_vector_free(F);
	gsl_matrix_free(sigtune);
	
	return 0;
}


