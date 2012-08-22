#include <stdlib.h>

void sort(double *x, int *idx, int n);
int EquipartitionYAxis(double *Dy, int n, int y, int *Qm);
int GetSuperclumpsPartition(double *Dx, int n, int *Qm, int k, int *Pm);
void ApproxOptimizeXAxis(double *Dx, double *Dy, int n, int *Qm, int q, 
			 int *Pm, int p, int x, double *I);
