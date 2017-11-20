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
#define erf_N1  0.80104398333355385655425471455314706359364542054
#define erro_GL     1.144827473
#define erro_GLr    4.204133489

/*O valor que aceitaremos como verdadeiro é erf(0.9083010)=0.80104398333355385655425471455314706359364542054.*/

/*Variáveis globais*/

/*Vetores de pesos e abcissas a serem usados na quadratura por Gauss-Legendre*/

float roots[15]={-0.987992518020485,-0.937273392400706,-0.848206583410427,-0.724417731360170,-0.570972172608539,-0.394151347077563,-0.201194093997434,0,0.201194093997434,0.394151347077563,0.570972172608539,0.724417731360170,0.848206583410427,0.937273392400706,0.987992518020485};
float weights[15]={0.0307532419961173,0.0703660474881081,0.1071592204671720,0.1395706779261540,0.1662692058169940,0.1861610000155620,0.1984314853271120,0.2025782419255610,0.1984314853271120,0.1861610000155620,0.1662692058169940,0.1395706779261540,0.1071592204671720,0.0703660474881081,0.0307532419961173};

/*Vetores de pesos e abcissas a serem usados na quadratura por Gauss-Kronrod*/

float abcissas_G_7[7]={0.949107912342759,0.741531185599394,0.405845151377397,0,-0.405845151377397,-0.741531185599394,-0.949107912342759};
float abcissas_K_15[15]={0.991455371120813,0.949107912342759,0.864864423359769,0.741531185599394,0.586087235467691,0.405845151377397,0.207784955007898,0,-0.207784955007898,-0.405845151377397,-0.586087235467691,-0.741531185599394,-0.864864423359769,-0.949107912342759,-0.991455371120813};
float pesos_G_7[7]={0.129484966168870,0.279705391489277,0.381830050505119,0.417959183673469,0.381830050505119,0.279705391489277,0.129484966168870};
float pesos_K_15[15]={0.022935322010529,0.063092092629979,0.104790010322250,0.140653259715525,0.169004726639267,0.190350578064785,0.204432940075298,0.209482141084728,0.204432940075298,0.190350578064785,0.169004726639267,0.140653259715525,0.104790010322250,0.063092092629979,0.022935322010529};

/*Protótipos de funções usadas ao longo do programa*/

float trapezio(float x1, float x2, int p);
void simpson(float x1, float x3, int p);
float gauss_legendre(int flag);
void gauss_kronrod(int flag);

void raiz(float x[15]);

float legendre_15(float x);
void teste();

/*Programa principal*/

int main(){
	int i=0;
	float j;

    /*teste();*/

    printf("Valor aceito\nerf(%9.7f)=%20.18f\n\n", N1, erf_N1);

    printf("\nMetodo do trapezio:\n");
	for(i=0; i<21; i++) trapezio(0, N1, i);

	printf("\nMetodo de Simpson:\n");
	for(i=0; i<21; i++) simpson(0, N1, i);

	printf("\nMetodo da quadratura de Gauss-Legendre:\n");
	j=gauss_legendre(FALSE);
    printf("erf(%9.7f)=%20.18f\n", N1, j);
    printf("Erro absoluto: %20.18f\n", fabs(j-erf_N1));
    printf("Erro de truncamento: -%fE-31\n", erro_GL); //o cálculo dessa constante se encontra nos comentários ao fim do código

	printf("\nMetodo da quadratura de Gauss-Kronrod:\n");
	gauss_kronrod(FALSE);

	printf("\nMetodo da quadratura de Gauss-Legendre (refinado):\n");
	j=gauss_legendre(TRUE);
    printf("erf(%9.7f)=%20.18f\n", N1, j);
    printf("Erro absoluto: %20.18f\n", fabs(j-erf_N1));
    printf("Erro de truncamento: -%fE-50\n", erro_GLr); //o cálculo dessa constante se encontra nos comentários ao fim do código

	printf("\nMetodo da quadratura de Gauss-Kronrod (refinado):\n");
    gauss_kronrod(TRUE);


	return 0;
}

float trapezio(float x1, float x2, int p){
	int i=0;
	float I, L;
	L=(x2-x1)/(pow(2, p+1)-1);  //em um segmento com n pontos há n-1 segmentos intermediários
	I=(exp(-x1*x1)+exp(-x2*x2))/2;
	for(i=1; i<pow(2, p+1)-1; i++) I += exp(-pow((x1+i*L), 2));
	I=2*L*I/sqrt(M_PI);
	printf("%d  %20.18f %30.28f %30.28f\n", p, I, fabs(I-erf_N1), pow(2,-(2*p+2)));
}

void simpson(float x1, float x3, int p){
	int i=0;
	float I, E, L, A, B, C;
	L=(x3-x1)/(2*(pow(3, p)));
    I=exp(-x1*x1)+exp(-x3*x3);

	for(i=1, A=0; i<=pow(3, p); i++)
        A += exp(-(x1+2*i*L)*(x1+2*i*L));
    A=2*A;

    for(i=1, B=0, E=0; i<=pow(3, p); i++){
        C=exp(-(x1+L+2*(i-1)*L)*(x1+L+2*(i-1)*L));
        E += C*(4*pow(x1+2*(i-1)*L, 4)-12*pow(x1+2*(i-1)*L, 2)+3); //derivada quarta de exp(-x²) no ponto à esquerda
        B += C;
    }
    B=4*B;

    I=I+A+B;
	I=((2*L*I)/(3*sqrt(M_PI)));
	E=fabs(8*E*pow(L, 5)/(90*sqrt(M_PI)));

	/*a expressão para o erro está em 25.4.5 do Abramowitz & Stegun;
	o ponto à esquerda é escolhido porque a derivada quarta é monotônica
	decrescente em todo intervalo e esse é um erro estimado por cima*/

	//printf("%d   %20.18f %30.28f\n", p, I, E);
	if(p<16) printf("%d  %20.18f    %30.28f     %30.28f\n", p, I, fabs(I-erf_N1), E);
	else{
        I=2*pow(3, p-16)*I;
        printf("%d  %20.18f    %30.28f     %30.28f\n", p, I, fabs(I-erf_N1), E);
	}
}

float gauss_legendre(int flag){
	int i,j;
	float I;
	if(flag==FALSE) raiz(roots);	//refina as raízes previamente escolhidas
    if(flag==TRUE){
        for(j=0, I=0; j<15; j++)
            for(i=0; i<4; i++) I += weights[j]*exp(-pow(((N1/8)*(roots[j]+2*i+1)), 2));
        return (N1*I)/(4*sqrt(M_PI));
    }
    else{
        for(i=0, I=0; i<15; i++) I += (weights[i])*exp(-pow((N1/2)*(roots[i]+1), 2));
        return (N1*I/sqrt(M_PI));
    }
}

void gauss_kronrod(int flag){
	int i, j;
	float I, A;
	I=A=0;
    float q[2]={0,0};
    float x[4]={0,0,0,0};
    float y[4]={0,0,0,0};

	if(flag==TRUE){

        for(i=0, I=0; i<4; i++){
            for(j=0; j<7; j++){
                I += pesos_G_7[j]*exp(-pow(((N1/8)*(abcissas_G_7[j]+2*i+1)), 2));
            }
            x[i]=(N1*I)/(4*sqrt(M_PI));             //cada um dos quatro G_7's
        }
        q[0]=(N1*I)/(4*sqrt(M_PI));
        for(i=0, I=0; i<4; i++){
            for(j=0; j<15; j++){
                I += pesos_K_15[j]*exp(-pow(((N1/8)*(abcissas_K_15[j]+2*i+1)), 2));
            }
            y[i]=(N1*I)/(4*sqrt(M_PI));             //cada um dos quatro K_15's
        }
        q[1]=(N1*I)/(4*sqrt(M_PI));
        for(i=0, A=0; i<4; i++){
            A += pow(fabs(200*(x[i]-y[i])), 1.5);   //soma dos erros de cada uma das quatro integrais
        }
    }

    else{
        for(i=0, I=0; i<7; i++) I += (pesos_G_7[i])*exp(-pow((N1/2)*(abcissas_G_7[i]+1), 2));
        q[0]=(N1*I/sqrt(M_PI));
        for(i=0, I=0; i<15; i++) I += (pesos_K_15[i])*exp(-pow((N1/2)*(abcissas_K_15[i]+1), 2));
        q[1]=(N1*I/sqrt(M_PI));
        A=pow(fabs(200*(q[0]-q[1])), 1.5);
    }
	printf("G_7: %40.38f\nK_15: %40.38lf\n", q[0], q[1]);
	printf("Erro absoluto: %20.18f\n", fabs(q[1]-erf_N1));
    printf("Erro de truncamento: %30.28f\n", A);
}

/*Refina as raizes do polinômio de Legendre
de ordem 15 estimadas a partir de valores
iniciais fornecidos por tabelas de raízes
do próprio polinômio*/

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

float legendre_15(float x){
	return (1/2048)*(9694845*pow(x, 15)-35102025*pow(x, 13)+50702925*pow(x, 11)-37182145*pow(x, 9)+14549535*pow(x, 7)-2909907*pow(x, 5)+255255*pow(x, 3)-6435*x);
}

void teste(){
	raiz(roots);
	int i;
	for(i=0; i<15; i++){
		printf("legendre_15(%20.18lf)=%20.18lf\n", roots[i], legendre_15(roots[i]));
	}
}

/*Alguns links utilizados:

http://www.wolframalpha.com/widgets/view.jsp?id=6c23cf54347d8e2dda52f41a918efb43
http://keisan.casio.com/exec/system/1180573449
https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights
https://pomax.github.io/bezierinfo/legendre-gauss.html*/

/*Algumas constantes que precisamos calcular no Wolfram Alpha:

Cálculo do erro do método de Gauss-Legendre (sem refinamento):
d^30(exp(-pow(((0.9083010/2)*(x+1)), 2)))/dx^30 (x=-1)=-1.05476*1E10=A => R_15 = [2³¹(15!)⁴/31(30!)³]*A = 1,144827473E-31;

Erro total de truncamento do método de Gauss-Legendre
calculado separadamente da rotina devido à precisão:
E=

Cálculo do erro do método de Gauss-Legendre (refinado):

d^30(exp(-pow(((0.9083010/8)*(x+2*0+1)), 2)))/dx^30 (x=-1)=-9.1486*1E-9=A  => R_15 = [2³¹(15!)⁴/31(30!)³]*A = −9,92981211E-50;
d^30(exp(-pow(((0.9083010/8)*(x+2*1+1)), 2)))/dx^30 (x=-1)=1.79427*1E-9=A => R_15 = 1,947485295E-50;
d^30(exp(-pow(((0.9083010/8)*(x+2*2+1)), 2)))/dx^30 (x=-1)=7.5953*1E-9=A => R_15 = 8,243873589E-50;
d^30(exp(-pow(((0.9083010/8)*(x+2*3+1)), 2)))/dx^30 (x=-1)=-4.11435*1E-9=A => R_15 = -4,465680263E-50;
[2³¹(15!)⁴/31(30!)³]=2,025660298×10⁵⁶/1,866294709×10⁹⁷=1,085391438E-41

Erro total de truncamento do método de Gauss-Legendre
calculado separadamente da rotina devido à precisão:
E=(−9,92981211+1,947485295+8,243873589−4,465680263)*1E-50=−4,204133489E-50*/
