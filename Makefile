CC = gcc
MPICC = mpicc
FLG = -O4
NAME = kmeansTest

all: kmeansTest.o kmeans.o cluster.h

	$(MPICC) $(FLG) kmeansTest.o kmeans.o -lm -o $(NAME)

kmeansTest.o: kmeansTest.c kmeans.h

	$(MPICC) $(FLG) kmeansTest.c -c -lm

kmeans.o: kmeans.c kmeans.h cluster.h

	$(MPICC) $(FLG) kmeans.c -c

clean:
	rm -f *.o *.out *.exe
	rm -f *.bin  

