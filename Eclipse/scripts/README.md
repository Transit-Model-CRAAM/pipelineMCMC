# 🧵 Ativando Multithread com OpenMP

Este projeto permite o uso de **OpenMP** para aceleração de cálculos paralelos via scripts em C. Abaixo estão os passos para compilar e testar o suporte a **multithread** no Windows, Linux ou macOS.

---

## 🪟 Windows

### 📁 1. Acesse a pasta `scripts`

```bash
cd Eclipse/scripts
```

### 🛠️ 2. Crie o script `compile_windows.bat`

<details>
<summary>Clique para ver o conteúdo</summary>

```bat
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
```
</details>

### 🚀 3. Execute o script

```bash
.\compile_windows.bat
```

### ✅ 4. Verifique multithread

```bash
.\multithreadCheck.exe
```

Resultado esperado:

```
OpenMP version: 201511
Max threads: 12
Thread 0 of 12 is running
Thread 1 of 12 is running
...
Thread 11 of 12 is running
```

---

## 🐧 Linux / 🍎 macOS

### 📁 1. Acesse a pasta `scripts`

```bash
cd Eclipse/scripts
```

### 🛠️ 2. Crie o script `compile_unix.sh`

<details>
<summary>Clique para ver o conteúdo</summary>

```bash
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
```
</details>

### 🚀 3. Torne o script executável e rode

```bash
chmod +x compile_unix.sh
./compile_unix.sh
```

### ✅ 4. Verifique a saída

Exemplo:

```
OpenMP version: 201511
Max threads: 8
Thread 0 of 8 is running
Thread 1 of 8 is running
...
```

---

## ⚠️ Requisitos

- GCC com suporte a OpenMP:
  - Windows: [MSYS2](https://www.msys2.org/) ou MinGW-w64
  - Linux/macOS: geralmente já incluso
- Arquivos `func.c` e `multithreadCheck.c` presentes em `Eclipse/scripts`

---

Pronto! Agora o ECLIPSE está pronto para rodar com paralelização e ganho de performance 🚀
