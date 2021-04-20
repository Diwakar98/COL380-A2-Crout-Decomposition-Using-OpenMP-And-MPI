all:
	gcc -fopenmp -O0 -o run main.c

run:
	./run 5 filename.txt 16 0