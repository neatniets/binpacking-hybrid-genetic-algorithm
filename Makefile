CC = gcc -g -pthread -O2
PF_PATH = parallel-foreach

main.out: main.o bin-packing.o chromosome.o bp-solution.o parallel-foreach.o
	$(CC) -o main.out main.o bin-packing.o chromosome.o bp-solution.o \
		parallel-foreach.o

main.o: main.c bin-packing.h bp-solution.h
	$(CC) -c main.c

parallel-foreach.o: $(PF_PATH).c $(PF_PATH).h
	$(CC) -c $(PF_PATH).c

bin-packing.o: bin-packing.c bin-packing.h bp-solution.h chromosome.h \
	$(PF_PATH).h
	$(CC) -c bin-packing.c

chromosome-test.out: chromosome.o chromosome-test.o bp-solution.o
	$(CC) -o chromosome-test.out chromosome.o chromosome-test.o \
		bp-solution.o

chromosome-test.o: chromosome-test.c chromosome.h
	$(CC) -c chromosome-test.c

chromosome.o: chromosome.c chromosome.h bp-solution.h
	$(CC) -c chromosome.c

bp-solution-test.out: bp-solution-test.o bp-solution.o
	$(CC) -o bp-solution-test.out bp-solution-test.o bp-solution.o

bp-solution-test.o: bp-solution-test.c bp-solution.h
	$(CC) -c bp-solution-test.c

bp-solution.o: bp-solution.h bp-solution.c
	$(CC) -c bp-solution.c

.PHONY : clean
clean:
	-rm *.o
	-rm *.out
	-rm *.gch
