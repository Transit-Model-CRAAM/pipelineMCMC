#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include "func.h"

int* criaEstrela(int lin, int col, int tamanhoMatriz, float raio, float intensidadeMaxima, float coeficienteHum, float coeficienteDois){
	int i, j;
	int *estrela = (int*) malloc (lin * col * sizeof(int));
	int index;
	float centro = tamanhoMatriz/2.0f;
	float raioSquared = raio * raio;

#pragma omp parallel
{
#pragma omp for collapse(2)
	for(i=0;i<lin;i++){
		for(j=0;j<col;j++){
			index = i*(col) + j;  // Fixed index calculation
			estrela[index] = 0;	
		}
	}
	float distanciaCentroSquared;
	float cosTheta;
	float dx, dy;
	float oneMinusCosTheta, oneMinusCosThetaSquared;

#pragma omp for collapse(2)
	for(j=0;j<col;j++){
		for(i=0;i<lin;i++){
			dx = i - centro;
			dy = j - centro;
			distanciaCentroSquared = dx*dx + dy*dy;

			if(distanciaCentroSquared <= raioSquared){
				cosTheta = sqrt(1.0f - distanciaCentroSquared/raioSquared);
				oneMinusCosTheta = 1.0f - cosTheta;
				oneMinusCosThetaSquared = oneMinusCosTheta * oneMinusCosTheta;

				index = i*(col) + j;  // Fixed index calculation
				estrela[index] = (int) (intensidadeMaxima * (1.0f - coeficienteHum * oneMinusCosTheta - coeficienteDois * oneMinusCosThetaSquared));
			}
		}
	}
}

	return estrela;
}

double curvaLuz(double x0, double y0, int tamanhoMatriz, double raioPlanetaPixel, double *estrelaManchada, double *kk, double maxCurvaLuz){
	double valor = 0;
	int i;
	double raioPlanetaPixelSquared = raioPlanetaPixel * raioPlanetaPixel;
	double dy, dx, distSquared;
	int row, col;
	double tamanhoMatrizDouble = (double)tamanhoMatriz;

#pragma omp parallel for reduction(+:valor) private(dy, dx, distSquared, row, col)
	for(i=0;i<tamanhoMatriz*tamanhoMatriz;i++){
		// Calculate row and column more efficiently
		row = (int)(kk[i] / tamanhoMatrizDouble);
		col = (int)(kk[i] - tamanhoMatrizDouble * row);

		// Calculate distance squared
		dy = row - y0;
		dx = col - x0;
		distSquared = dy*dy + dx*dx;

		if(distSquared > raioPlanetaPixelSquared){
			valor += estrelaManchada[i];
		}
	}

	valor = valor/maxCurvaLuz;
	return valor;
}

double curvaLuzLua(double x0, double y0, double xm, double ym, double rMoon, int tamanhoMatriz, double raioPlanetaPixel, double *estrelaManchada, double *kk, double maxCurvaLuz){
	double valor = 0;
	int i;
	double raioPlanetaPixelSquared = raioPlanetaPixel * raioPlanetaPixel;
	double rMoonSquared = rMoon * rMoon;
	double dy1, dx1, distSquared1;  // For planet
	double dy2, dx2, distSquared2;  // For moon
	int row, col;
	double tamanhoMatrizDouble = (double)tamanhoMatriz;

#pragma omp parallel for reduction(+:valor) private(dy1, dx1, distSquared1, dy2, dx2, distSquared2, row, col)
	for(i=0;i<tamanhoMatriz*tamanhoMatriz;i++){
		// Calculate row and column more efficiently
		row = (int)(kk[i] / tamanhoMatrizDouble);
		col = (int)(kk[i] - tamanhoMatrizDouble * row);

		// Calculate distance squared for planet
		dy1 = row - y0;
		dx1 = col - x0;
		distSquared1 = dy1*dy1 + dx1*dx1;

		// Calculate distance squared for moon
		dy2 = row - ym;
		dx2 = col - xm;
		distSquared2 = dy2*dy2 + dx2*dx2;

		// Check if point is outside both planet and moon
		if(distSquared1 > raioPlanetaPixelSquared && distSquared2 > rMoonSquared){
			valor += estrelaManchada[i];
		}
	}

	// Normalizacao
	valor = valor/maxCurvaLuz;
	return valor;
}

// Function to free memory allocated by criaEstrela
void liberaEstrela(int *estrela) {
    if (estrela != NULL) {
        free(estrela);
    }
}