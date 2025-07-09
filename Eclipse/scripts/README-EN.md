# 🧵 Enabling Multithreading with OpenMP

This project uses **OpenMP** to accelerate calculations with C parallel scripts. Below are the instructions to compile and test **multithread support** on **Windows**, **Linux**, or **macOS**.

---

## 🪟 Windows

### 📁 1. Navigate to the `scripts` folder

```bash
cd Eclipse/scripts
```

### 🛠️ 2. Create a file named `compile_windows.bat`

<details>
<summary>Click to expand the content</summary>

```bat
@echo off
REM === Compile with OpenMP support ===

REM Compile 32-bit version
gcc -fopenmp -Wall -std=c11 -c -fPIC func.c -o func32.o -O3
gcc -fopenmp -shared func32.o -o func32.dll -O3

REM Compile 64-bit version
gcc -fopenmp -Wall -std=c11 -c -fPIC func.c -o func64.o -O3
gcc -fopenmp -shared func64.o -o func64.dll -O3

REM Compile multithread check program
gcc -fopenmp -Wall -std=c11 multithreadCheck.c -o multithreadCheck.exe

echo.
echo === Compilation complete! ===
echo Running multithread check:
multithreadCheck.exe
```
</details>

### 🚀 3. Run the script

```bash
.\compile_windows.bat
```

### ✅ 4. Run multithread check

```bash
.\multithreadCheck.exe
```

Expected output:

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

### 📁 1. Navigate to the `scripts` folder

```bash
cd Eclipse/scripts
```

### 🛠️ 2. Create a file named `compile_unix.sh`

<details>
<summary>Click to expand the content</summary>

```bash
#!/bin/bash

echo "Compiling with OpenMP for Linux/macOS..."

# Compile object file with OpenMP support
gcc -fopenmp -Wall -std=c11 -fPIC -c func.c -o func.o -O3

# Create shared library
gcc -fopenmp -shared func.o -o libfunc.so -O3

# Compile multithread check program
gcc -fopenmp -Wall -std=c11 multithreadCheck.c -o multithreadCheck -O3

echo "Compilation complete!"
echo "Running multithread check:"
./multithreadCheck
```
</details>

### 🚀 3. Make the script executable and run it

```bash
chmod +x compile_unix.sh
./compile_unix.sh
```

### ✅ 4. Verify output

Example:

```
OpenMP version: 201511
Max threads: 8
Thread 0 of 8 is running
Thread 1 of 8 is running
...
```

---

## ⚠️ Requirements

- GCC with OpenMP support:
  - **Windows**: [MSYS2](https://www.msys2.org/) or MinGW-w64
  - **Linux/macOS**: Usually preinstalled
- `func.c` and `multithreadCheck.c` files must be in `Eclipse/scripts`

---

You're all set! Your project is now configured to run with multithreading and performance gains 🚀
