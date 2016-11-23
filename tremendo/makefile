CC=g++
LIBDIR=/usr/local/lib
CFLAGS=-I.:/usr/local/include -L$(LIBDIR) -lgsl -lgslcblas
DEPS = test.h structs.h definiciones.h tremendo.h
OBJ= capture.o knockout.o onetrans.o tremendo.o
OBJtest = test.o testsub.o 
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

tremendo: $(OBJ) 
	$(CC) -o tremendo $(OBJ) $(CFLAGS) 


testexe: $(OBJtest) 
	$(CC) -o testexe $(OBJtest) $(CFLAGS) 