compiler = g++
flags = -std=c++20 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
flagsForOpenMP = -L/opt/homebrew/opt/libomp/lib -lomp

main: main.cpp solver.h
	$(compiler) $(flags) -o main main.cpp $(flagsForOpenMP)

test: main
	./main
