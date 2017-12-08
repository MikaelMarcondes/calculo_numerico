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

double C_S=332; //velocidade do som, em m/s
double A=M_PI*(pow(0.15, 2)/4); //área da seção transversa, em m2
//double A=0; //área da seção transversa, em m2, particula puntiforme
double rho=1.2; //densidade em kg/m3
double M[9]={0, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3};
double C_D[9]={1.6, 1.8, 2.1, 4.8, 5.8, 5.4, 4.8, 4.2, 3.9};
double g=9.79; //gravidade em m/s2
double V0=700; //velocidade inicial em m/s
double m=50; //massa em kg

/*Protótipos de funções usadas ao longo do programa*/
double lagr_pol(double a, double b, double c, double A, double B, double C, double x);
double drag_coefficient(double s);
double speed(double v[2]);
void runge_kutta_fehlberg(double v[2], int flaf, double dt);
double f(double v_1, int flag, double v_2);

/*Programa principal*/

int main(){
    double x[2]={0,0};

    double a;
    printf("Entre com o ângulo inicial (inteiro de 1 a 6):");
    scanf("%lf",&a);
    a *= 10*M_PI*N2/180;
    double v[2]={V0*cos(a), V0*sin(a)};
    double z[2]={0,0}; //vetor para guardar temporariamente as velocidades antigas

    int i;
    double dt=1E-3;

    FILE *fp;
    fp=fopen("trajetoria.txt", "w");

    if(fp == NULL)
        exit(-1);
    for(i=0; x[1]>=0; i++){
        z[0]=x[0];
        z[1]=x[1];
        runge_kutta_fehlberg(v, TRUE, dt); //com flag=TRUE altera a primeira componente da velocidade
        runge_kutta_fehlberg(v, FALSE, dt); //com flag=FALSE altera a segunda componente da velocidade
        x[0] += dt*v[0];
        x[1] += dt*v[1];

        fprintf(fp, "%6.4lf   %6.4lf\n", z[0], z[1]);
    }
    /*fprintf(fp, "\nAlcance do projetil: %6.4lf m\n", z[0]);
    fprintf(fp, "Tempo de voo: %6.4lf s\n", i*dt);
    fprintf(fp, "Velocidade de impacto: %6.4lf m/s\n", speed(v));*/
    fclose(fp);
}

double lagr_pol(double a, double b, double c, double A, double B, double C, double x){
    return ((A*(x-b)*(x-c)/((a-b)*(a-c)))+(B*(x-a)*(x-c)/((b-a)*(b-c)))+(C*(x-a)*(x-b)/((c-a)*(c-b))));
}


double drag_coefficient(double s){
    int i,j,k;
    i=j=k=0;
    s=s/C_S;    //numero de Mach
    for(i=0; i<9; i=i+3){
        j=i+1;
        k=i+2;
        if(M[i]<=s && s<=M[k]){
            return lagr_pol(M[i], M[j], M[k], C_D[i], C_D[j], C_D[k], s);
        }
    }
}

double speed(double v[2]){
    return (pow((pow(v[0], 2)+pow(v[1], 2)), 0.5));
}

void runge_kutta_fehlberg(double v[2], int flag, double dt){
    int i=0;
    if(flag==TRUE) i=0;
    else i=1;

    double k[6];

    k[1]=dt*f(v[i], flag, v[((i+1)%2)]);
    k[2]=dt*f(v[i]+(k[1]/4), flag, v[((i+1)%2)]);
    k[3]=dt*f(v[i]+(3*k[1]/32)+(9*k[2]/32), flag, v[((i+1)%2)]);
    k[4]=dt*f(v[i]+(1932*k[1]/2197)-(7200*k[2]/2197)+(7296*k[3]/2197), flag, v[((i+1)%2)]);
    k[5]=dt*f(v[i]+(439*k[1]/216)-8*k[2]+(3680*k[3]/513)-(845*k[4]/4104), flag, v[((i+1)%2)]);
    k[6]=dt*f(v[i]-(8*k[1]/27)+2*k[2]-(3544*k[3]/2565)-(1859*k[4]/4104)-(11*k[5]/40), flag, v[((i+1)%2)]);
    v[i] += (16*k[1]/135)+(6656*k[3]/12825)+(28561*k[4]/56430)-(9*k[5]/50)+(2*k[6]/55);
}

double f(double v_1, int flag, double v_2){
    double s;
    double V[2]={v_1, v_2};
    s=speed(V);
    if(flag==TRUE) return (-A*rho*drag_coefficient(s)*v_1*s/(2*m));
    else return (-A*rho*drag_coefficient(s)*v_1*s/(2*m)-g);
}
