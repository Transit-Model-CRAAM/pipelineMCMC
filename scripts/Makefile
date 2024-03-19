CC=clang
FLAGS=-g -Wall -std=c++17 -arch arm64 -lomp
all: func64.dylib
func64.dylib: func64.o
	$(CC) $(FLAGS) -shared -arch x86_64 func64.o -o func64.dylib -O3
func64.o: func.c
	clang -std=c11 -arch x86_64 -Wall -Wextra -pedantic -c -fPIC func.c -o func64.o -O3