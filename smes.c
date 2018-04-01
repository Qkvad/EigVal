/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: System of masses with elastic strings
 *
 *
 * BUILD: gcc smes.c -o smes_exe -lblas -llapack -lm
 * RUN  : ./smes_exe
 *
 ====================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

int main(){
	integer i, n=4, N=n*n;
	doublereal M[16]={2.0, 0.0, 0.0, 0.0,
					  0.0, 5.0, 0.0, 0.0, 
					  0.0, 0.0, 3.0, 0.0, 
					  0.0, 0.0, 0.0, 6.0};
	doublereal A[16]={24.0, -9.0, -5.0, 0.0, 
					  -9.0, 22.0, -8.0, -5.0, 
					  -5.0, -8.0, 25.0, -7.0, 
					  0.0, -5.0, -7.0, 18.0};
	doublereal* work=malloc(N*sizeof(doublereal));
	doublereal* D=malloc(n*sizeof(doublereal));
	doublereal* u=malloc(n*sizeof(doublereal));
	doublereal* p=malloc(n*sizeof(doublereal));
	
	for(i=0;i<n;i++) M[i+i*n]=1.0/sqrt(M[i+i*n]);
	
	char trans='N'; 
	doublereal alpha=1.0,beta=0.0;
	
	dgemm_(&trans,&trans,&n,&n,&n,&alpha,M,&n,A,&n,&beta,work,&n);
	dgemm_(&trans,&trans,&n,&n,&n,&alpha,work,&n,M,&n,&beta,A,&n);

	char jobz='V',uplo='L';
	integer info, ldw=3*n-1;
	dsyev_(&jobz,&uplo,&n,A,&n,D,work,&ldw,&info);

	//for(i=0;i<N;i++) printf("%lf\n",A[i]);
	//printf("svojstvene vrijednosti\n");
	//for(i=0;i<n;i++) printf("%lf\n",D[i]);

	//spremi u work umnozak M*A
	dgemm_(&trans,&trans,&n,&n,&n,&alpha,M,&n,A,&n,&beta,work,&n);

	integer j,k,inc=1;
	printf("\n");
	for(i=0;i<n;i++){
		printf("x%ld=(",i+1);
		for(k=0;k<n;k++){
				p[k]=A[k+i*n];
			}
		dgemv_(&trans,&n,&n,&alpha,M,&n,p,&inc,&beta,u,&inc);
		for(j=0;j<n;j++){
			if(j==n-1) printf("%lf",u[j]);
			else printf("%lf, ",u[j]);
		}
		printf(")*exp(%lf *t*i)\n",sqrt(D[i]));
	}
	printf("\n");
	
	free(work);free(D);free(u);free(p);
return 0;
}
