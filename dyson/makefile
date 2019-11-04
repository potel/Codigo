#
#  makefile 
#
#echo $LDFLAGS
CC=g++

DIR := build
DIR_RH7 := build/rh7

################
#  C++ files   #
################
SRC := dyson.cpp lagrange.cpp spectralFunction.cpp volume.cpp hartreeFock.cpp pot.cpp read_parameters.cpp RedSelfEnergy.cpp kMesh.cpp potPara.cpp imaginaryForm.cpp expInt.cpp spinOrbit.cpp asy.cpp gaussInteg.cpp gauss16.cpp disperse.cpp twoFermi.cpp surfaceTF.cpp surVolume.cpp standard.cpp surfaceSTD.cpp io.cpp eigen.cpp density.cpp numerical.cpp surfaceGeneral.cpp sphericalB.cpp legendre.cpp

OBJ_RH7 := $(SRC:%.cpp=$(DIR_RH7)/%.o)
OBJ := $(SRC:%.cpp=$(DIR)/%.o)

##################################################################################################################
#These includes and libraries are machine-dependent, make sure they correspond to yours if you want this to work##
##################################################################################################################

# Includes
#LOCAL_INCLUDE = -I./usr/local/include/
LOCAL_INCLUDE = -I./
GSL_CFLAGS = $(shell /usr/bin/gsl-config --cflags)
#BNB_INCLUDE = -I/opt/boost_numeric_bindings-v1/
BNB_INCLUDE = -I/usr/local/numeric_bindings-v1
BNB_INCLUDE_COMP = -I/usr/include/boost-numeric-bindings/
BOOST = -I/usr/local/boost_1_63_0
EIGEN_INCLUDE = -I/usr/local/eigen3
# Libraries
LIBDIR=-L/usr/local/lib
LAPACK_LIBS = -llapack
ATLAS_LIBS = -L/usr/local/Atlas/lib 
ATLAS_LIBS_RH7 = -L/usr/lib64/atlas -lsatlas 
ATLAS_LIBS_COMP = -L/usr/local/Atlas/lib  -larmadillo -latlas	 -pg
GSL_LIBS = $(shell /usr/bin/gsl-config --libs)

CFLAGS= -w -c -O3 $(BNB_INCLUDE) $(LOCAL_INCLUDE) $(GSL_CFLAGS) $(EIGEN_INCLUDE) $(BOOST) -fpermissive
LINKOPTION_RH7 = $(ATLAS_LIBS_RH7) $(LAPACK_LIBS) $(GSL_LIBS)
#LINKOPTION =  -llapack $(GSL_LIBS) $(LIBDIR) $(ATLAS_LIBS_COMP) $(LAPACK_LIBS) $(ATLAS_LIBS)
LINKOPTION =  -I.:/usr/local/include $(LIBDIR) -llapack -lgsl -lgslcblas -larmadillo -pg
#Note that the option: -std=c++03 is a weird option I had to add to compile the calli code on my laptop. It may or may not be needed for general machines.
#I had to put the option in because my compiler did not like the various uses of "getline" in calli.cpp for some reason...
CFLAGS_COMP= -w -c -O3 $(BNB_INCLUDE) $(LOCAL_INCLUDE) $(GSL_CFLAGS) -fpermissive -std=c++03

##########################################################################################
#programs
##########################################################################################

green: $(DIR)/green.o $(OBJ) 
	$(CC) -o green $(DIR)/green.o $(OBJ) $(LINKOPTION)

calli: $(DIR)/calli.o $(OBJ)
	$(CC) -o calli $(DIR)/calli.o $(OBJ) $(LINKOPTION)

green_rh7: $(DIR_RH7)/green_rh7.o $(OBJ_RH7) $(OBJ_RH7_FOR) 
	$(CC) -o green_rh7 $(DIR_RH7)/green_rh7.o $(OBJ_RH7_FOR) $(OBJ_RH7)  $(LINKOPTION_RH7)

###########################################################################################
###########################################################################################


$(DIR_RH7)/%.o: %.cpp | $(DIR_RH7)
	$(CC) $(CFLAGS) $< -o $@

$(DIR)/%.o: %.cpp | $(DIR)
	$(CC) $(CFLAGS_COMP) $< -o $@

$(DIR_RH7)/%.o : %.f | $(DIR_RH7)
	$(FC) $(FFLAGS) $< -o $@

$(DIR)/%.o : %.f | $(DIR)
	$(FC) $(FFLAGS) $< -o $@

$(DIR)/green.o: green.cpp | $(DIR)
	$(CC) $(CFLAGS_COMP) $< -o $@

$(DIR_RH7)/green_rh7.o: green.cpp | $(DIR_RH7)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(DIR)/*.o green
cleanrh7:
	rm $(DIR_RH7)/*.o green_rh7
cleanall:
	rm $(DIR)/*.o $(DIR_RH7)/*.o green_rh7 green

