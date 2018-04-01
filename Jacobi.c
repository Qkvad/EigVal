/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: Jacobi method for spectral decomposition - Ris matrix
 *
 *
 * BUILD: gcc Jacobi.c -o Jacobi_exe -lblas -llapack -lm
 * RUN  : ./Jacobi_exe
 *
 ====================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

doublereal norma(integer n, doublereal* A){
	int i,j;
	doublereal sum=0.0;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if(i!=j) sum+=A[i+j*n]*A[i+j*n];
		}
	}
	return sqrt(sum);
}

int sgn(doublereal x){
	if(x>=0) return 1;
	else return -1;
}	

void jacobi_sd(integer n, doublereal* A, doublereal tol){
	integer p,q,k;
	doublereal tau,t,s,c,app,apq,aqq,pom,norm;
	doublereal* work=malloc(n*sizeof(doublereal));	
	char nrm='F';
	//tol*=dlange_(&nrm,&n,&n,A,&n,work);
	norm=norma(n,A); 
	while(norm>tol*dlange_(&nrm,&n,&n,A,&n,work)){
		for(p=0;p<n-1;p++){
			for(q=p+1;q<n;q++){
				if(A[p+q*n]!=0){
					tau=(A[q+q*n]-A[p+p*n])/(2*A[p+q*n]);
					t=sgn(tau)/(fabs(tau)+sqrt(1+pow(tau,2)));
					c=1/(sqrt(1+pow(t,2)));
					s=t*c;
				}
				else{
					c=1.0; s=0.0;	
				}
				app=A[p+p*n];
				apq=A[p+q*n];
				aqq=A[q+q*n];
				app=app-t*apq;
				aqq=aqq+t*apq;
				for(k=0;k<n;k++){
					pom=A[k+p*n];
					A[k+p*n]=c*pom-s*A[k+q*n];
					A[k+q*n]=s*pom+c*A[k+q*n];
					A[p+k*n]=A[k+p*n];
					A[q+k*n]=A[k+q*n];			
				}
				A[p+q*n]=0.0; A[q+p*n]=0.0;
				A[p+p*n]=app; A[q+q*n]=aqq;
			}	
		}
		norm=norma(n,A); printf("%lf\n",norm);	
	}
}

int main(){
	integer n=10,N=n*n,i,j;
	doublereal tol=1e-15;
	doublereal* A=malloc(N*sizeof(doublereal));	

	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			A[i+j*n]=1.0/(16.0-2.0*i-2.0*j+3.0);
		}
	}

	printf("\n norme:\n");
	jacobi_sd(n,A,tol);
	printf("\n svojstvene vrijednosti:\n");
	for(i=0;i<n;i++) printf("%lf \n", A[i+i*n]);

	free(A);
return 0;
}
