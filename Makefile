all:
	make 2
	make 3
	make 4
	make 5
	make 6
	make gs
	make jac


clean:
	rm -rvf *~
	rm -rvf solved?
	rm -rvf gs-omp jacobi-omp
	clear

push:
	git push https://github.com/yairdaon/HPC15-HW3.git

# jacobi with omp
jacobi-omp:jacobi-omp.c
	gcc -fopenmp -Wall jacobi-omp.c -lm -o jacobi-omp
jac: jacobi-omp
	./jacobi-omp 500


# GS with OMP
gs-omp: gs-omp.c
	gcc -fopenmp -Wall gs-omp.c -lm -o gs-omp	
gs: gs-omp
	./gs-omp 500





# the bugs!!!
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

solved6: omp_solved6.c
	gcc -fopenmp omp_solved6.c -o solved6
6: solved6
	./solved6

