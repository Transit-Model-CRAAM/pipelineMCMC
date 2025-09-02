#!/bin/bash

echo "Compilando com OpenMP para Linux/macOS..."

# Compilar objeto com suporte a OpenMP
gcc -fopenmp -Wall -std=c11 -fPIC -c func.c -o func.o -O3

# Criar biblioteca compartilhada
gcc -fopenmp -shared func.o -o libfunc.so -O3

# Compilar verificador
gcc -fopenmp -Wall -std=c11 multithreadCheck.c -o multithreadCheck -O3

echo "Compilação concluída!"
echo "Executando verificação de multithread:"
./multithreadCheck