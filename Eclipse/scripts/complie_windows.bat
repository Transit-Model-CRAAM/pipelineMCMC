@echo off
REM === Compilação com suporte a OpenMP ===

REM Compilar versão 32 bits
gcc -fopenmp -Wall -std=c11 -c -fPIC func.c -o func32.o -O3
gcc -fopenmp -shared func32.o -o func32.dll -O3

REM Compilar versão 64 bits
gcc -fopenmp -Wall -std=c11 -c -fPIC func.c -o func64.o -O3
gcc -fopenmp -shared func64.o -o func64.dll -O3

REM Compilar programa de verificação de multithread
gcc -fopenmp -Wall -std=c11 multithreadCheck.c -o multithreadCheck.exe

echo.
echo === Compilação concluída! ===
echo Executando verificação de multithread:
multithreadCheck.exe