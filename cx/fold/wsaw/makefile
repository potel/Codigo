FCOMP = gfortran -c
FLINK = gfortran

OBJ = wsaw.o

wsaw.o:		wsaw.for
	$(FCOMP) wsaw.for
wsaw: wsaw.o 
	$(FLINK) $(OBJ) -o wsaw

