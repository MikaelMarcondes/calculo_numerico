/*     Nome: Mikael Jorge Marcondes de Oliveira - Nº USP: 9302351      */
/*           Nome: Rafael Bicudo Ribeiro - Nº USP: 9083010	       */
/*   Disciplina: MAP 0214 - Cálculo Numérico com Aplicações em Física  */
/* Instituto de Física da Universidade de São Paulo - Segundo Semestre */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Constantes recorrentes*/

#define	N1	0.9083010
#define	N2	0.9302351
#define FALSE	0
#define TRUE	1

/*O valor que aceitaremos como verdadeiro é erf(0.9083010)=0.80104398333355385655425471455314706359364542054.*/

/*Variáveis globais*/

/*Vetores de pesos e abcissas a serem usados na quadratura por Gauss-Legendre*/

double roots[15]={-0.987992518020485,-0.937273392400706,-0.848206583410427,-0.724417731360170,-0.570972172608539,-0.394151347077563,-0.201194093997434,0,0.201194093997434,0.394151347077563,0.570972172608539,0.724417731360170,0.848206583410427,0.937273392400706,0.987992518020485};
double weights[15]={0.0307532419961173,0.0703660474881081,0.1071592204671720,0.1395706779261540,0.1662692058169940,0.1861610000155620,0.1984314853271120,0.2025782419255610,0.1984314853271120,0.1861610000155620,0.1662692058169940,0.1395706779261540,0.1071592204671720,0.0703660474881081,0.0307532419961173};

/*Vetores de pesos e abcissas a serem usados na quadratura por Gauss-Kronrod*/

double abcissas_G_7[7]={0.949107912342759,0.741531185599394,0.405845151377397,0,-0.405845151377397,-0.741531185599394,-0.949107912342759};
double abcissas_K_15[15]={0.991455371120813,0.949107912342759,0.864864423359769,0.741531185599394,0.586087235467691,0.405845151377397,0.207784955007898,0,-0.207784955007898,-0.405845151377397,-0.586087235467691,-0.741531185599394,-0.864864423359769,-0.949107912342759,-0.991455371120813};
double pesos_G_7[7]={0.129484966168870,0.279705391489277,0.381830050505119,0.417959183673469,0.381830050505119,0.279705391489277,0.129484966168870};
double pesos_K_15[15]={0.022935322010529,0.063092092629979,0.104790010322250,0.140653259715525,0.169004726639267,0.190350578064785,0.204432940075298,0.209482141084728,0.204432940075298,0.190350578064785,0.169004726639267,0.140653259715525,0.104790010322250,0.063092092629979,0.022935322010529};

/*Protótipos de funções usadas ao longo do programa*/

double trapezio(double x1, double x2, int p);
double simpson(double x1, double x3, int p);
double gauss_legendre(int flag, double a, double b, double alpha, double beta);
double gauss_kronrod(int flag, double a, double b, double alpha, double beta, double q[2]);

double f(double x);
void raiz(double x[15]);

double legendre_15(double x);
void teste();

double transf_linear(double x, double a, double b, double alpha, double beta);

/*double d_legendre_15(double x);
void pesos();/*
/*Programa principal*/

int main(){
	int i=0;
	double j, k;
	double a[2]={0,0};
	
	printf("\nMetodo do trapezio:\n");
	for(i=0; i<21; i++) printf("%d	%20.18lf \n", i, trapezio(0, N1, i));

	printf("\nMetodo de Simpson:\n");
	for(i=0; i<7; i++) printf("%d	%20.18lf \n", i, simpson(0, N1, i));

	printf("\nMetodo da quadratura de Gauss-Legendre:\n");
	printf("erf(%lf)=%20.18lf\n", N1, gauss_legendre(FALSE, 0, 0, 0, 0)); //os valores depois de FALSE são irrelevantes

	printf("\nMetodo da quadratura de Gauss-Kronrod:\n");
	gauss_kronrod(FALSE, 0, 0, 0, 0, a);
	printf("G_7: %20.18lf\nK_15: %20.18lf\n", a[0], a[1]);

	printf("\nMetodo da quadratura de Gauss-Legendre (refinado):\n");
	for(i=0; i<4; i++) j += gauss_legendre(TRUE,((-1)+(i/2)),((-1/2)+(i/2)), -1, 1);
	printf("erf(%lf)=%20.18lf\n", N1, j);

	printf("\nMetodo da quadratura de Gauss-Kronrod (refinado):\n");
	for(i=0, j=0, k=0, a[0]=0, a[1]=0; i<4; i++){
		gauss_kronrod(TRUE,((-1)+(i/2)),((-1/2)+(i/2)), -1, 1, a);
		j += a[0];
		k += a[1];
	}
	printf("G_7: %20.18lf\nK_15: %20.18lf\n", j, k);

	return 0;
}

double trapezio(double x1, double x2, int p){
	int i=0;
	double I, L;
	I=0;
	L=(x2-x1)/pow(2, p+1);
	for(i=0, x2=x1+L; i<pow(2, p+1); i++){
		I += (f(x1)+f(x2))*L/2;
		x1 += L;
		x2 += L;
	}
	return I;
}

double simpson(double x1, double x3, int p){
	int i=0;
	double I, L, x2;
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

double gauss_legendre(int flag, double a, double b, double alpha, double beta){
	int i=0;
	double I=0;
	double t=0;
	if(flag==FALSE) raiz(roots);	//refina as raízes previamente escolhidas
	for(i=0; i<15; i++){
		if(flag==TRUE) t=transf_linear(roots[i], a, b, alpha, beta);
		else t=roots[i];
		I += (weights[i])*exp(-pow((N1/2)*(t+1), 2));
	}
	if(flag==TRUE) return (((b-a)/(beta-alpha))*(N1*I/sqrt(M_PI)));
	return (N1*I/sqrt(M_PI));
}

double gauss_kronrod(int flag, double a, double b, double alpha, double beta, double q[2]){
	int i=0;
	double I, t;
	I=t=0;

	for(i=0; i<7; i++){
		if(flag==TRUE) t=transf_linear(abcissas_G_7[i], a, b, alpha, beta);
		else t=abcissas_G_7[i];
		I += (pesos_G_7[i])*exp(-pow((N1/2)*(t+1), 2));
	}
	if(flag==TRUE) q[0]=(((b-a)/(beta-alpha))*(N1*I/sqrt(M_PI)));
	else q[0]=(N1*I/sqrt(M_PI));

	for(i=0, I=0; i<15; i++){
		if(flag==TRUE) t=transf_linear(abcissas_K_15[i], a, b, alpha, beta);
		else t=abcissas_K_15[i];
		I += (pesos_K_15[i])*exp(-pow((N1/2)*(t+1), 2));
	}
	if(flag==TRUE) q[1]=(((b-a)/(beta-alpha))*(N1*I/sqrt(M_PI))); //(b-a)/(beta-alpha)=1/4
	else q[1]=(N1*I/sqrt(M_PI));

	//printf("Erro estimado: %20.18lf \n", pow(200*fabs(g_7-k_15), 1.5));
}

double f(double x){
	return ((2/sqrt(M_PI)) * exp(-x*x));
}

/*Refina as raizes do polinômio de Legendre
de ordem 15 estimadas a partir de valores
iniciais fornecidos por tabelas de raízes
do próprio polinômio*/

void raiz(double x[15]){
	int i, j;
	double y, a, b;
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

double legendre_15(double x){
	return (1/2048)*(9694845*pow(x, 15)-35102025*pow(x, 13)+50702925*pow(x, 11)-37182145*pow(x, 9)+14549535*pow(x, 7)-2909907*pow(x, 5)+255255*pow(x, 3)-6435*x);
}

void teste(){
	raiz(roots);
	int i;
	for(i=0; i<15; i++){
		printf("legendre_15(%20.18lf)=%20.18lf\n", roots[i], legendre_15(roots[i]));
	}
}

double transf_linear(double x, double a, double b, double alpha, double beta){
	return ((((b-a)*x)+(a*beta-b*alpha))/(beta-alpha));
}

/*double d_legendre_15(double x){
	return (1/2048)*((15*9694845)*pow(x, 14)-(13*35102025)*pow(x, 12)+(11*50702925)*pow(x, 10)-(9*37182145)*pow(x, 8)+(7*14549535)*pow(x, 6)-(5*2909907)*pow(x, 4)+(3*255255)*pow(x, 2)-6435);
}
void pesos(){
	int i;
	for(i=0; i<15; i++){
		weights[i]=2/((1-pow(roots[i], 2))*pow(d_legendre_15(roots[i]), 2));
	}
}*/
/*Alguns links utilizados:

http://www.wolframalpha.com/widgets/view.jsp?id=6c23cf54347d8e2dda52f41a918efb43
http://keisan.casio.com/exec/system/1180573449
https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights
https://pomax.github.io/bezierinfo/legendre-gauss.html*/
