/**
 * @file func.h
 * @brief Header file for Eclipse simulation functions
 * 
 * This file contains declarations for functions used in the Eclipse simulation,
 * including star creation, light curve calculation, and memory management.
 */

#ifndef FUNC_H
#define FUNC_H

/**
 * @brief Creates a star with limb darkening
 * 
 * @param lin Number of rows in the star matrix
 * @param col Number of columns in the star matrix
 * @param tamanhoMatriz Size of the matrix (typically equal to lin and col)
 * @param raio Radius of the star in pixels
 * @param intensidadeMaxima Maximum intensity at the center of the star
 * @param coeficienteHum First limb darkening coefficient
 * @param coeficienteDois Second limb darkening coefficient
 * @return int* Pointer to the star matrix (must be freed with liberaEstrela)
 */
int* criaEstrela(int lin, int col, int tamanhoMatriz, float raio, float intensidadeMaxima, float coeficienteHum, float coeficienteDois);

/**
 * @brief Calculates the light curve for a planet transiting a star
 * 
 * @param x0 X-coordinate of the planet
 * @param y0 Y-coordinate of the planet
 * @param tamanhoMatriz Size of the matrix
 * @param raioPlanetaPixel Radius of the planet in pixels
 * @param estrelaManchada Pointer to the star matrix
 * @param kk Pointer to the pixel indices array
 * @param maxCurvaLuz Maximum value of the light curve (for normalization)
 * @return double Normalized flux value
 */
double curvaLuz(double x0, double y0, int tamanhoMatriz, double raioPlanetaPixel, double *estrelaManchada, double *kk, double maxCurvaLuz);

/**
 * @brief Calculates the light curve for a planet with a moon transiting a star
 * 
 * @param x0 X-coordinate of the planet
 * @param y0 Y-coordinate of the planet
 * @param xm X-coordinate of the moon
 * @param ym Y-coordinate of the moon
 * @param rMoon Radius of the moon in pixels
 * @param tamanhoMatriz Size of the matrix
 * @param raioPlanetaPixel Radius of the planet in pixels
 * @param estrelaManchada Pointer to the star matrix
 * @param kk Pointer to the pixel indices array
 * @param maxCurvaLuz Maximum value of the light curve (for normalization)
 * @return double Normalized flux value
 */
double curvaLuzLua(double x0, double y0, double xm, double ym, double rMoon, int tamanhoMatriz, double raioPlanetaPixel, double *estrelaManchada, double *kk, double maxCurvaLuz);

/**
 * @brief Frees memory allocated by criaEstrela
 * 
 * @param estrela Pointer to the star matrix to be freed
 */
void liberaEstrela(int *estrela);

#endif /* FUNC_H */