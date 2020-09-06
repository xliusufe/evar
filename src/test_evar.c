#include <math.h> 			
#include <stdio.h>  		
#include <stdlib.h> 		
#include <Rinternals.h>  	
#include <float.h>  		
#include <string.h> 		
#include "rbs_hdh.h"


double EstEvar1(double *y, double beta0, double tau, const int n, int maxstep, double eps){
	int i,step = 0;
	double tmp, beta=0.0, w, yw, tau1 = 1.0-tau;

	while (step<maxstep){
		step++;
		tmp = yw = 0.0;
		for(i=0;i<n;i++){
			w	= (y[i]>beta0)?tau:tau1;
			tmp	+= w;
			yw	+= y[i]*w;
		}
		if(tmp<DBL_EPSILON) return 0.0;
		beta	= yw/tmp;
		if(fabs(beta-beta0)<eps) break;
		else	beta0 = beta;
	}
	return beta;
}

void EstEvar(double *beta, double *beta0, double *x, double *y, double tau, const int n, int p, int maxstep, double eps){
	int i,j,step = 0;
	double tmp,bnorm, tausqrt=sqrt(tau), tausqrt1=sqrt(1-tau);
	double *w, *xw, *yw, *Q, *R, *qy, *invR;
	xw		= (double*)malloc(sizeof(double)*n*p);
	yw		= (double*)malloc(sizeof(double)*n);
	w		= (double*)malloc(sizeof(double)*n);
	Q 		= (double*)malloc(sizeof(double)*n*p);  
	R 		= (double*)malloc(sizeof(double)*p*p);  
	invR 	= (double*)malloc(sizeof(double)*p*p);	
	qy 		= (double*)malloc(sizeof(double)*p);   	

	while (step<maxstep){
		step++;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++)
				tmp	+= x[j*n+i]*beta0[j];
			w[i]	= (y[i]>tmp)?tausqrt:tausqrt1;
		}
		for(j=0;j<p;j++)
			for(i=0;i<n;i++)
				xw[j*n+i]	= x[j*n+i]*w[i]; 
		for(i=0;i<n;i++)
			yw[i]	= y[i]*w[i]; 

		QRDecompN(Q, R, xw, n, p);
		LowTriangularInv(invR, p, R);
		for(j=0;j<p;j++){
			tmp = 0.0; 
			for(i=0;i<n;i++) tmp += Q[j*n+i]*yw[i];
			qy[j] = tmp;
		}
		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp = 0.0;
			for(i=j;i<p;i++)	tmp	+= invR[i*p+j]*qy[i];
			beta[j]	= tmp;
			bnorm	+= (tmp-beta0[j])*(tmp-beta0[j]);
		}
		if(sqrt(bnorm)<eps) break;
		else{
			for(j=0;j<p;j++)	beta0[j] = beta[j];
		}
	}
	free(xw);
	free(yw);
	free(w);
	free(Q);
	free(R);
	free(invR);
	free(qy);	
}

double Test_Evar1(const double *x, const double *y, double beta, double tau, const int n, int p){
	int i,j,k;
	double *ytau,*epshat;
	double tmp, trSigma2=0.0, Tn=0.0;

	ytau  	= (double*)malloc(sizeof(double)*n);
	epshat  = (double*)malloc(sizeof(double)*n);

	for(i=0;i<n;i++){
		tmp = y[i] - beta;
		ytau[i] = (tmp>0.0)?tau: (1-tau);
		epshat[i]	= tmp;
	}			
	for(i=1;i<n;i++){
		for(j=0;j<i;j++){			
			tmp = 0.0;
			for(k=0;k<p;k++)	tmp	+= x[k*n+i]*x[k*n+j];
			tmp		*= epshat[i]*epshat[j]*ytau[i]*ytau[j];
			trSigma2+= tmp*tmp;
			Tn	 	+= tmp; 
		}				
	}
	Tn 			*= 2.0/n/(n-1);
	trSigma2 	*= 2.0/n/(n-1);	
	Tn			*= n/sqrt(2.0*trSigma2);

	free(ytau);	
	free(epshat);
	return Tn;
}

double Test_Evar(const double *x1, const double *x, const double *y, double *beta, double tau, const int n, int p, int q){
	int i,j,k;
	double *ytau,*epshat;
	double tmp, trSigma2=0.0, Tn=0.0;

	ytau  	= (double*)malloc(sizeof(double)*n);
	epshat  = (double*)malloc(sizeof(double)*n);

	for(i=0;i<n;i++){
		tmp = y[i];
		for(j=0;j<q;j++)
			tmp	-= x1[j*n+i]*beta[j];
		ytau[i] = (tmp>0.0)?tau: (1-tau);
		epshat[i]	= tmp;
	}			
	for(i=1;i<n;i++){
		for(j=0;j<i;j++){			
			tmp = 0.0;
			for(k=0;k<p;k++)	tmp	+= x[k*n+i]*x[k*n+j];
			tmp		*= epshat[i]*epshat[j]*ytau[i]*ytau[j];
			trSigma2+= tmp*tmp;
			Tn	 	+= tmp; 
		}				
	}
	Tn 			*= 2.0/n/(n-1);
	trSigma2 	*= 2.0/n/(n-1);	
	Tn			*= n/sqrt(2.0*trSigma2);

	free(ytau);	
	free(epshat);
	return Tn;
}

SEXP EST_EVAR_(SEXP X_, SEXP Y_, SEXP BETA0_, SEXP PARA_INT_, SEXP PARA_DOUBLE_)
{
	// dimensions
	int *para 	= INTEGER(PARA_INT_);
	int n     	= para[0];
	int p     	= para[1];
	int maxstep	= para[2];

	double *para1 	= REAL(PARA_DOUBLE_);
	double tau   	= para1[0];
	double eps    	= para1[1];

	// Pointers
	double *y		= REAL(Y_);
	double *x		= REAL(X_);
	double *beta0 	= REAL(BETA0_);

	// Outcome
	SEXP _output, _beta, _r_names;
	PROTECT(_output 	= allocVector(VECSXP, 	1));
	PROTECT(_r_names   	= allocVector(STRSXP, 	1));
	PROTECT(_beta 		= allocVector(REALSXP, 	p));
	

	if(p==1)
		REAL(_beta)[0]	= EstEvar1(y, beta0[0], tau, n, maxstep, eps);
	else			
		EstEvar(REAL(_beta), beta0, x, y, tau, n, p, maxstep, eps);
	

	SET_STRING_ELT(_r_names, 	0,	mkChar("beta"));
	SET_VECTOR_ELT(_output,		0, 	_beta);
	setAttrib(_output, 			R_NamesSymbol, _r_names); 

	UNPROTECT(3);
	return _output;
}

SEXP TEST_EVAR_(SEXP X1_, SEXP X2_, SEXP Y_, SEXP BETA0_, SEXP PARA_INT_, SEXP PARA_DOUBLE_)
{
	// Pointers
	int *para 	= INTEGER(PARA_INT_);
	int n     	= para[0];
	int p     	= para[1];
	int q       = para[2];
	int maxstep	= para[3];

	double *para1 	= REAL(PARA_DOUBLE_);
	double tau   	= para1[0];
	double eps    	= para1[1];

	// Pointers
	double *y		= REAL(Y_);
	double *x1		= REAL(X1_);	
	double *x2		= REAL(X2_);
	double *beta0 	= REAL(BETA0_);

	// Outcome
	SEXP _output, _beta, _Tn, _r_names;
  	PROTECT(_output		= allocVector(VECSXP, 	2));
  	PROTECT(_r_names   	= allocVector(STRSXP, 	2));
	PROTECT(_Tn 		= allocVector(REALSXP, 	1));
	PROTECT(_beta 		= allocVector(REALSXP, 	q));
  	

	if(q==1){
		REAL(_beta)[0]	= EstEvar1(y, beta0[0], tau, n, maxstep, eps);
		REAL(_Tn)[0]	= Test_Evar1(x2, y, REAL(_beta)[0], tau, n, p);
	}
	else{			
		EstEvar(REAL(_beta), beta0, x1, y, tau, n, q, maxstep, eps);
		REAL(_Tn)[0]	= Test_Evar(x1, x2, y, REAL(_beta), tau, n, p, q);
	}

  	

  	SET_STRING_ELT(_r_names,	0,	mkChar("Tn"));
	SET_STRING_ELT(_r_names,	1,	mkChar("beta"));
	SET_VECTOR_ELT(_output, 	0, 	_Tn);
	SET_VECTOR_ELT(_output, 	1, 	_beta);
	setAttrib(_output, 			R_NamesSymbol,	_r_names); 

	UNPROTECT(4);
	return _output;
}






