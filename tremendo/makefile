CC=g++
OBJ= tremendo.o capture.o knockout.o onetrans.o raditrans.o racah_algebra.o
LIBDIR=/usr/local/lib
#CFLAGS= -std=c++11    -I.:/usr/local/include -L$(LIBDIR) -lgsl -lgslcblas -larmadillo
#CFLAGS= -std=c++14     -lgsl -lgslcblas -larmadillo
#CFLAGS= -std=c++11 -lgsl -lgslcblas -larmadillo
CFLAGS= -std=c++11    -I.:/opt/homebrew/Cellar/armadillo/12.6.4/include -L$(LIBDIR) -L//opt/homebrew/Cellar/armadillo/12.6.4/lib -lgsl -lgslcblas -larmadillo
DEPS= structs.h definiciones.h tremendo.h
tremendo: $(OBJ) 
	$(CC)  -o $@  $(OBJ) $(CFLAGS)



%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o

