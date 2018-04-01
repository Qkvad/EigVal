/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: main angles and vectors - SVD
 *
 *
 * BUILD: gcc ang.c -o ang_exe -lblas -llapack -lm
 * RUN  : ./ang_exe
 *
 ====================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

void ispisMatrice(doublereal *A,integer m,integer n){
	int i,j;
	printf("\n");
	for(i=0;i<m;i++){
		for(j=0;j<n;j++)
			printf("%f ",A[i+j*m]);
		printf("\n");
	}
	printf("\n");
}

int main(){
	int i;
	integer m=4,n=2,N=n*n;
	doublereal A[]={1.0,0.0,0.0,0.0,
					1.0,1.0,1.0,1.0};
	doublereal B[]={1.0,-1.0,1.0,-1.0,
					0.0,1.0,0.0,1.0};
	
	//QRf od A i B
	doublereal *TAUA=malloc(n*sizeof(doublereal));
	doublereal *TAUB=malloc(n*sizeof(doublereal));
	doublereal *work=malloc(N*sizeof(doublereal));
	integer info;
	dgeqrf_(&m,&n,A,&m,TAUA,work,&N,&info);
	dgeqrf_(&m,&n,B,&m,TAUB,work,&N,&info);

	//A=Q iz QRf od A
	dorgqr_(&m,&n,&n,A,&m,TAUA,work,&N,&info);
	//B=Q iz QRf od B
	dorgqr_(&m,&n,&n,B,&m,TAUB,work,&N,&info);
	
	//C=A^t *B
	doublereal C[N];
	char transA='T',transB='N';
	doublereal alpha=1.0,beta=0.0;
	dgemm_(&transA,&transB,&n,&n,&m,&alpha,A,&m,B,&m,&beta,C,&n);

	//SVD dio
	char jobu='A',jobvt='A';
	integer ldw=5*n;
	doublereal S[n],U[N],VT[N],work1[ldw];
	dgesvd_(&jobu,&jobvt,&n,&n,C,&n,S,U,&n,VT,&n,work1,&ldw,&info);

	printf("\ncos kutova = ");
	for(i=0;i<n;i++) printf("%lf   ",S[i]);
	printf("\n");

	//X=A*U    Y=B*V
	doublereal X[m*n],Y[m*n];
	dgemm_(&transB,&transB,&m,&n,&n,&alpha,A,&m,U,&n,&beta,X,&m);	
	dgemm_(&transB,&transA,&m,&n,&n,&alpha,B,&m,VT,&n,&beta,Y,&m);
	
	printf("\nVektori x");
	ispisMatrice(X,m,n);
	printf("\nVektori y");
	ispisMatrice(Y,m,n);

	//min arccos(ddot(x,y)/||x||2 * ||y||2)
	doublereal x[m],y[m];	
	integer inc=1,k=0;
	for(i=m;i<m*n;i++){
		x[k]=X[i];
		y[k]=Y[i];
		k++;
	}
	doublereal br=ddot_(&m,x,&inc,y,&inc);
	doublereal norma_x=dnrm2_(&m,x,&inc);
	doublereal norma_y=dnrm2_(&m,y,&inc);
	doublereal kut=acos(br/(norma_x*norma_y));
	printf("kut izmedu potprostora = %lf\n",kut);

	free(TAUA);free(TAUB);free(work);
return 0;
}
