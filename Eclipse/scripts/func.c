#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

int* criaEstrela(int lin, int col, int tamanhoMatriz, float raio, float intensidadeMaxima, float coeficienteHum, float coeficienteDois){
	int i, j;
	int *estrela = (int*) malloc (lin * col * sizeof(int*));
	int index;

#pragma omp parallel
{
#pragma omp for collapse(2)
	for(i=0;i<lin;i++){
		for(j=0;j<col;j++){
			index = i*(lin) + j;
			estrela[index] = 0;	
		}
	}
	float distanciaCentro;
	float cosTheta;
	
#pragma omp for collapse(2)
	for(j=0;j<col;j++){
		for(i=0;i<lin;i++){
			distanciaCentro = sqrt(pow(i-tamanhoMatriz/2,2) + pow(j-tamanhoMatriz/2,2));
			if(distanciaCentro <= raio){
				cosTheta = sqrt(1-pow(distanciaCentro/raio,2));
				index = i*(lin) + j;
				estrela[index] = (int) (intensidadeMaxima * (1 - coeficienteHum * (1 - cosTheta) - coeficienteDois * (pow(1 - cosTheta,2))));
			}
		}
	}
}

	return estrela;
}

double curvaLuz(double x0, double y0, int tamanhoMatriz, double raioPlanetaPixel, double *estrelaManchada, double *kk, double maxCurvaLuz){
	double valor = 0;
	int i;
	
#pragma omp parallel for reduction(+:valor)
	for(i=0;i<tamanhoMatriz*tamanhoMatriz;i++){
		if(pow((kk[i]/tamanhoMatriz-y0),2) + pow((kk[i]-tamanhoMatriz*floor(kk[i]/tamanhoMatriz)-x0),2) > pow(raioPlanetaPixel,2)){
			valor += estrelaManchada[i];
		}
	}
	
	valor = valor/maxCurvaLuz;

//((kk/tamanhoMatriz-y0)**2+(kk-tamanhoMatriz*np.fix(kk/tamanhoMatriz)-x0)**2 <= raioPlanetaPixel**2)	
//(self.estrelaManchada*plan,dtype=float)/maxCurvaLuz
	
	return valor;
}

double curvaLuzLua(double x0, double y0, double xm, double ym, double rMoon, int tamanhoMatriz, double raioPlanetaPixel, double *estrelaManchada, double *kk, double maxCurvaLuz){
	double valor = 0;
	int i;
	
#pragma omp parallel for reduction(+:valor)
	for(i=0;i<tamanhoMatriz*tamanhoMatriz;i++){
		if((pow((kk[i]/tamanhoMatriz-y0),2) + pow((kk[i]-tamanhoMatriz*floor(kk[i]/tamanhoMatriz)-x0),2) > pow(raioPlanetaPixel,2)) && (pow((kk[i]/tamanhoMatriz-ym),2) + pow((kk[i]-tamanhoMatriz*floor(kk[i]/tamanhoMatriz)-xm),2) > pow(rMoon,2))){
			valor += estrelaManchada[i];
		}
	}
	
	valor = valor/maxCurvaLuz;

//((kk/tamanhoMatriz-y0)**2+(kk-tamanhoMatriz*np.fix(kk/tamanhoMatriz)-x0)**2 <= raioPlanetaPixel**2)	
//(self.estrelaManchada*plan,dtype=float)/maxCurvaLuz
	
	return valor;
}
