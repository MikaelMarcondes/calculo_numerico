/*     Nome: Mikael Jorge Marcondes de Oliveira - Nº USP: 9302351      */
/*           Nome: Rafael Bicudo Ribeiro - Nº USP: 9083010	       */
/*   Disciplina: MAP 0214 - Cálculo Numérico com Aplicações em Física  */
/* Instituto de Física da Universidade de São Paulo - Segundo Semestre */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Constantes definidas para serem usadas*/

#define	N1	0.9083010
#define	N2	0.9302351

/*Variáveis globais*/

float roots[15]={-0.987992518020485,-0.937273392400706,-0.848206583410427,-0.724417731360170,-0.570972172608539,-0.394151347077563,-0.201194093997434,0,0.201194093997434,0.394151347077563,0.570972172608539,0.724417731360170,0.848206583410427,0.937273392400706,0.987992518020485};	//zeros do polinômio de Legendre de ordem 15

float weights[15]={0.0307532419961173,0.0703660474881081,0.1071592204671720,0.1395706779261540,0.1662692058169940,0.1861610000155620,0.1984314853271120,0.2025782419255610,0.1984314853271120,0.1861610000155620,0.1662692058169940,0.1395706779261540,0.1071592204671720,0.0703660474881081,0.0307532419961173};

/*Protótipos de funções usadas ao longo do programa*/

float trapezio(float x1, float x2, int p);
float simpson(float x1, float x3, int p);
float gauss_quad(float x[15]);
float f(float x);
float d_legendre_15(float x);
float legendre_15(float x);
void raiz(float x[15]);
void teste();
void pesos();

/*Programa principal*/

int main(){
	int i=0;
	raiz(roots);
	/*for(i=0; i<21; i++)
		printf("%d	%20.18lf \n", i, simpson(0, N1, i));*/
	printf("%20.18lf", gauss_quad(roots));
	return 0;
}

float trapezio(float x1, float x2, int p){
	int i=0;
	float I, L;
	I=0;
	L=(x2-x1)/pow(2, p+1);
	for(i=0, x2=x1+L; i<pow(2, p+1); i++){
		I += (f(x1)+f(x2))*L/2;
		x1 += L;
		x2 += L;
	}
	return I;
}

float simpson(float x1, float x3, int p){
	int i=0;
	float I, L, x2;
	I=0;
	L=(x3-x1)/(2*(pow(3, p)));
	for(i=0, x2=x1+L, x3=x1+2*L; i<pow(3, p); i++){
		I += (f(x1)+4*f(x2)+f(x3))*(L/3);
		x1 += 2*L;
		x2 += 2*L;
		x3 += 2*L;
	}
	return I;
}

float gauss_quad(float x[15]){
	int i=0;
	float I=0;
	for(i=8; i<9; i++){
		I += (weights[i])*exp(-pow((N1/2)*(x[i]+1), 2));
	}
	return (N1*I/sqrt(M_PI));
}

float f(float x){
	return ((2/sqrt(M_PI)) * exp(-x*x));
}

float legendre_15(float x){
	return (1/2048)*(9694845*pow(x, 15)-35102025*pow(x, 13)+50702925*pow(x, 11)-37182145*pow(x, 9)+14549535*pow(x, 7)-2909907*pow(x, 5)+255255*pow(x, 3)-6435*x);
}

float d_legendre_15(float x){
	return (1/2048)*((15*9694845)*pow(x, 14)-(13*35102025)*pow(x, 12)+(11*50702925)*pow(x, 10)-(9*37182145)*pow(x, 8)+(7*14549535)*pow(x, 6)-(5*2909907)*pow(x, 4)+(3*255255)*pow(x, 2)-6435);
}

/*Encontra as raizes do polinômio de Legendre
de ordem 15 a partir de valores iniciais fornecidos
por tabelas de raízes do próprio polinômio*/

void raiz(float x[15]){
	int i, j;
	float y, a, b;
	for(i=0; i<7; i++){
		y=x[i];
		a=x[i]-0.005;
		b=x[i]+0.005;
		while(fabs(legendre_15(x[i]))>1E-6){
			if(legendre_15(a)*legendre_15(x[i])<0){
				b=y;
				y=(a+b)/2;
			}
			else{
				a=y;
				y=(a+b)/2;
			}
		}
		x[i]=y;
	}
	for(i=8; i<15; i++){
		j=14-i;
		x[i]=-x[j];
	}
}

/*Testa se os zeros da tabela
são uma boa representação*/

void teste(){
	raiz(roots);
	int i;
	for(i=0; i<15; i++){
		printf("legendre_15(%25.23lf)=%25.23lf\n", roots[i], legendre_15(roots[i]));
	}
}

void pesos(){
	int i;
	for(i=0; i<15; i++){
		//weights[i]=2/((1-pow(roots[i], 2))*pow(d_legendre_15(roots[i]), 2));
		printf("%25.23lf\n", pow(d_legendre_15(roots[i]), 2));
	}
}
