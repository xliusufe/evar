#ifndef RBS_H_INCLUDED
#define RBS_H_INCLUDED

void sortN(int *ind0, double *x, int n, int dd);

int LowTriangularInv(double *B, int n, double *A);

void QRDecompN(double *E, double *R, double *x, int n, int p);

#endif // RBS_H_INCLUDED
