#include <stdio.h>  	
#include <stdlib.h> 	
#include <string.h> 	
#include <math.h>
#include <float.h>
#include "rbs_hdh.h"

void sortN(int *ind0, double *x, int n, int dd){
	int i, j, MaxInd, d, *ind;
	double tmp;
	ind = (int*)malloc(sizeof(int)*n);
	for(i=0;i<n;i++) ind[i] = i;

	d = (dd==n?dd-1:dd);
	for(i=0;i<d;i++)
	{
		tmp = x[0]; MaxInd = ind[0];
		for(j=1;j<n-i;j++)
		{
			if(x[j]<tmp)
			{
				x[j-1] = x[j];
				x[j] = tmp;
				ind[j-1] = ind[j];
				ind[j] = MaxInd;
			}
			else
			{
				tmp = x[j];
				MaxInd = ind[j];
			}
		}	
	}
	for(j=0;j<dd;j++) ind0[j] = ind[n-j-1];
	free(ind);
}

int LowTriangularInv(double *B, int n, double *A){ 
	// Input:
	// A is a lower triangular matrix
	//
	// Output:
	// B = inv(A)
	//	
	int i,j,k;
	const double EPS=DBL_EPSILON;
	for(i=0;i<n;i++)
		if(fabs(A[i*n+i])<EPS)	return(0);
	for(i=0;i<n;i++)	B[i*n+i] = 1;
	for(j=1;j<n;j++)
		for(i=0;i<j;i++)	B[j*n+i] = 0;

	for(i=n-1;i>=0;i--)//rows
	{
		if(fabs(A[i*n+i]-1)>EPS)
			for(j=i;j<n;j++)
				B[j*n+i] = B[j*n+i]/A[i*n+i];
		if(i>0)
		{
			for(j=i;j<n;j++)// columns
				for(k=0;k<i;k++)// rows
					B[j*n+k] = B[j*n+k] - A[i*n+k]*B[j*n+i];
		}
	}
	return(1);
}

void QRDecompN(double *E, double *R, double *x, int n, int p){
	// Input:
	// X is a p*n matrix
	//
	// Output:
	// R is a p*p lower triangular matrix
	// E is a p*n matrix satisfying E*t(E) = I_p
	//		
	double *Z, *znorm;
	double  tmp, tmp1;
	int i,j, k;
	
	Z = (double*)malloc(sizeof(double)*n*p);
	znorm = (double*)malloc(sizeof(double)*p);

	// calculate the first column
	tmp = 0;
	for(i=0;i<n;i++){
		Z[i] = x[i];
		tmp += Z[i]*Z[i];		
	}
	znorm[0] = sqrt(tmp);
	tmp = 0;
	for(i=0;i<n;i++){
		E[i] = x[i]/znorm[0];
		tmp += E[i]*x[i];
	}
	R[0] = tmp;

	//iteration from j=1...p	
	for(j=1;j<p;j++){		
		for(k=0;k<j;k++){
			tmp=0;	for(i=0;i<n;i++) tmp += E[k*n+i]*x[j*n+i];
			R[j*p+k] = tmp;
		}
		tmp1 = 0;
		for(i=0;i<n;i++){
			tmp = 0; for(k=0;k<j;k++) tmp += R[j*p+k]*E[k*n+i];
			Z[j*n+i] = x[j*n+i] - tmp;
			tmp1 += pow(Z[j*n+i],2);
		}
		znorm[j] = sqrt(tmp1);
		tmp1 = 0;
		for(i=0;i<n;i++) E[j*n+i] = Z[j*n+i]/znorm[j];
		for(i=0;i<n;i++) tmp1 += E[j*n+i]*x[j*n+i];
		R[j*p+j] = tmp1;
	}
	free(Z); free(znorm);
}

