/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: Graph partitioning
 *
 *
 * BUILD: gcc gp.c -o gp_exe -lblas -llapack -lm
 * RUN  : ./gp_exe
 *
 ====================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

void ispis_matrice( doublereal *A, integer n )
{
	int i, j;
	printf( "\n\n" );
	for( i = 0; i < n; i++ )
	{
		for( j = 0; j < n; j++ )
			printf( " %f ", A[ i + j*n ] );
		printf( "\n" );
	}
	printf( "\n\n" );
}


int main(){
	integer n=7, N=n*n, i, j;
	doublereal L[]={0,2,3,4,0,0,0,
					2,0,0,7,1,0,0,
					3,0,0,3,0,2,1,
					4,7,3,0,0,0,0,
					0,1,0,0,0,7,3,
					0,0,2,0,7,0,5,
					0,0,1,0,3,5,0};
	doublereal* D=calloc(N,sizeof(doublereal));
	for(i=0;i<n;i++){		
		for(j=0;j<n;j++)
			if(i!=j) {
				L[i+i*n]+=L[j+i*n];
				if(L[j+i*n]!=0.0) L[j+i*n]*=-1;
			}
		D[i+i*n]=1/sqrt(L[i+i*n]);
	}
	//LN=D*L*D	
	doublereal LN[N],pom[N];
	char trans='N';
	doublereal alpha=1.0,beta=0.0;
	dgemm_(&trans,&trans,&n,&n,&n,&alpha,D,&n,L,&n,&beta,pom,&n);
	dgemm_(&trans,&trans,&n,&n,&n,&alpha,pom,&n,D,&n,&beta,LN,&n);

	char jobz='V',range='I',uplo='L',cmach='S';
	integer il=2,iu=2, m=il-iu+1,lwork=8*n,iwork[5*n],ifail[N],info;
	doublereal vl, vu, abstol=2*dlamch_(&cmach), w1[n],w2[n],z1[n*m],z2[n*m],work[lwork];
	
	dsyevx_(&jobz,&range,&uplo,&n,L,&n,&vl,&vu,&il,&iu,&abstol,&m,w1,z1,&n,work,&lwork,iwork,ifail,&info);
	dsyevx_(&jobz,&range,&uplo,&n,LN,&n,&vl,&vu,&il,&iu,&abstol,&m,w2,z2,&n,work,&lwork,iwork,ifail,&info);

	printf("\nu_2=(");
	for(i=0;i<n;i++) {
		if(i!=(n-1)) printf("%lf, ",z1[i]);
		else printf("%lf) za lambda_2=%lf\n",z1[i],w1[0]);
	}	
	int V_2[n];
	int br=0;
	printf("V_1 = { ");
	for(i=0;i<n;i++){
		if(z1[i]<0) {V_2[br]=i+1;br++;}
		else printf("%d ",i+1);
	}
	printf("}\nV_2 = { ");
	for(i=0;i<br;i++) printf("%d ",V_2[i]);
	printf("}\n");

	printf("\nu_N,2=(");
	for(i=0;i<n;i++) {
		if(i!=(n-1)) printf("%lf, ",z2[i]);
		else printf("%lf) za lambda_N,2=%lf\n",z2[i],w2[0]);
	}
	
	integer inc=1;
	doublereal y[n];
	dgemv_(&trans,&n,&n,&alpha,D,&n,z2,&inc,&beta,y,&inc);
	
	br=0;
	printf("V_N,1 = { ");
	for(i=0;i<n;i++){
		if(z2[i]<0) {V_2[br]=i+1;br++;}
		else printf("%d ",i+1);
	}
	printf("}\nV_N,2 = { ");
	for(i=0;i<br;i++) printf("%d ",V_2[i]);
	printf("}\n");	

	free(D);
return 0;
}
