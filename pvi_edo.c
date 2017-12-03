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

/*Variaveis globais*/

double C_S=332 //velocidade do som, em m/s
double A=M_PI*(pow(0.15, 2)/4); //área da seção transversa, em m2
double rho=1.2 //densidade em kg/m3
double M[9]={0, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3};
double C_D[9]={1.6, 1.8, 2.1, 4.8, 5.8, 5.4, 4.8, 4.2, 3.9};

/*Protótipos de funções usadas ao longo do programa*/

/*Programa principal*/

int main(){

}

double lagr_pol(double a, double b, double c, double A, double B, double C, double x){
    return ((A*(x-b)*(x-c)/((a-b)*(a-c)))+(B*(x-a)*(x-c)/((b-a)*(b-c)))+(C*(x-a)*(x-b)/((c-a)*(c-b))));
}

double drag_force(double v[2]){
    double s;
    double d;
    s=speed(v);
    d=drag_coefficient(s);
    return (A*d*rho*pow(s,2)/2);
}

double drag_coefficient(double s){
    int i,j,k;
    i=j=k=0;
    s=s/C_s;    //numero de Mach
    for(i=0; i<9; i=i+3){
        j=i+1;
        k=i+2
        if(M[i]<=s && s<=M[k]){
            return lagr_pol(M[i], M[j], M[k], C_D[i], C_D[j], C_D[k], s);
        }
    }
}

double speed(double v[2]){
    return (pow(pow(v[0], 2)+pow(v[1], 2)), 0.5);
}
