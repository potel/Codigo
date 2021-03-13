CC=g++
OBJ= tremendo.o capture.o knockout.o onetrans.o raditrans.o racah_algebra.o
#OBJ= tremendo_a_paco.o capture.o knockout.o onetrans.o raditrans.o racah_algebra.o
LIBDIR=/usr/local/lib
CFLAGS= -std=c++14 -I.:/usr/local/include -L$(LIBDIR) -lgsl -lgslcblas -larmadillo
DEPS= structs.h definiciones.h tremendo.h
tremendo: $(OBJ) 
	$(CC)  -o $@  $(OBJ) $(CFLAGS)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o

