/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: Procrustes problem with SVD
 *
 *
 * BUILD: gcc pro.c -o pro_exe -lblas -llapack -lm
 * RUN  : ./pro_exe
 *
 ====================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

void ispisMatrice(doublereal *A,integer n){
	int i,j;
	printf("\n");
	for(i=0;i<n;i++){
		for(j=0;j<n;j++)
			printf("%f ",A[i+j*n]);
		printf("\n");
	}
	printf("\n");
}

int main(){
	int i;
	integer m=4, n=2, N=n*n;
	doublereal A[]={1.2, 2.9, 5.2, 6.8, 
					2.1, 4.3, 6.1, 8.1};
	doublereal B[]={1.0, 3.0, 5.0, 7.0, 
					2.0, 4.0, 6.0, 8.0};
	doublereal C[N];
	
	char transA='T',transB='N';
	doublereal alpha=1.0,beta=0.0;
	dgemm_(&transA,&transB,&n,&n,&m,&alpha,A,&m,B,&m,&beta,C,&n);

	char jobu='A',jobvt='A';
	integer ldw=5*n, info;
	doublereal S[n],U[N],VT[N],work[ldw];
	dgesvd_(&jobu,&jobvt,&n,&n,C,&n,S,U,&n,VT,&n,work,&ldw,&info);

	printf("Singularne vrijednosti = {");
	for(i=0;i<n-1;i++) printf("%lf,",S[i]);
	printf("%lf}\n",S[n-1]);

	doublereal Q[N];
	dgemm_(&transB,&transB,&n,&n,&n,&alpha,U,&n,VT,&n,&beta,Q,&n);
	ispisMatrice(Q,n);

	beta=-1.0;
	dgemm_(&transB,&transB,&m,&n,&n,&alpha,A,&m,Q,&n,&beta,B,&m);
	char norm='F';
	doublereal norma=dlange_(&norm,&m,&n,B,&m,work);
	printf("||AQ-B||_f = %lf\n",norma);
	
return 0;
}
