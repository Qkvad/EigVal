/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: Least squares problem with svd
 *
 *
 * BUILD: gcc nksvd.c -o nksvd_exe -lblas -llapack -lm
 * RUN  : ./nksvd_exe
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
	integer n=2,m=10,N=m*n,c=n*n,a=m*m;
	doublereal A[N];
	doublereal b[]={3.5,4.9,6.8,9.3,10.9,13.4,15.1,16.7,19.0,21.2};
	for(i=0;i<m;i++) A[i]=1.0;
	for(i=m;i<2*m;i++) A[i]=i-9;
	//for(i=0;i<N;i++) printf("%lf \n",A[i]);

	char jobu='A',jobvt='A';
	integer ldw=5*m,info;
	doublereal work[ldw];
	doublereal S[n],U[a],VT[c];
	dgesvd_(&jobu,&jobvt,&m,&n,A,&m,S,U,&m,VT,&n,&work,&ldw,&info);
	doublereal* sigma=calloc(c,sizeof(doublereal));
	for(i=0;i<n;i++) sigma[i+i*n]=S[i];

	/*printf("\nsigma=\n");
	ispisMatrice(S,n,n);
	printf("\nU=\n");
	ispisMatrice(U,m,n);
	printf("\nV=\n");
	ispisMatrice(VT,n,n);*/
	
	//y=alpha * U^t * b + beta * y
	char trans='T';
	integer inc=1;
	doublereal alpha=1.0,beta=0.0;
	doublereal* y=malloc(n*sizeof(doublereal));
	dgemv_(&trans,&m,&n,&alpha,U,&m,b,&inc,&beta,y,&inc);
	
	//sigma = sigma ^ -1
	for(i=0;i<n;i++) sigma[i+i*n]=1.0/sigma[i+i*n];
	
	//C=V*sigma
	trans='N';
	doublereal* C=malloc(n*n*sizeof(doublereal));
	dgemm_(&trans,&trans,&n,&n,&n,&alpha,VT,&n,sigma,&n,&beta,C,&n);

	//a=C*y
	doublereal* aa=malloc(n*sizeof(doublereal));
	dgemv_(&trans,&n,&n,&alpha,C,&n,y,&inc,&beta,aa,&inc);

	printf("\n p(x)= %.4lf + %.4lfx \n",aa[0],aa[1]);

	free(sigma); free(y);free(C);
return 0;
}
