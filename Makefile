#COMPILER = gcc -fopenmp
#EXECUTABLES = omp_solved2 omp_solved3 omp_solved4 omp_solved5 omp_solved6

solved2: omp_solved2.c
	gcc  -fopenmp omp_solved2.c -o solved2

2: solved2
	./solved2



solved3: omp_solved3.c
	gcc -fopenmp omp_solved3.c -o solved3

3: solved3
	./solved3





solved4: omp_solved4.c
	gcc -fopenmp omp_solved4.c -o solved4

4: solved4
	./solved4



solved5: omp_solved5.c
	gcc -fopenmp omp_solved5.c -o solved5

5: solved5
	./solved5



solved6: solved6.c
	gcc -fopenmp omp_solved6.c -o solved6

6: solved
	./solved6



clean:
	rm -rvf *~
	rm -rvf solved?
	clear

push:
	git push https://github.com/yairdaon/HPC15-HW3.git