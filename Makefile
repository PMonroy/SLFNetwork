CC_MONO=g++ -Wall -mcmodel=medium

CC=$(CC_MONO)

LIBS= -lnetcdf_c++
RM=rm -rf

all: grid_construction.out lagrangian_engine.out matrix_construction.out

parameters.o: parameters.cpp ioutil.o parameters.h date.h
	$(CC) -c parameters.cpp

grid_construction.o: grid_construction.cpp parameters.o velocity.o ioutil.o 
	$(CC) -c grid_construction.cpp

lagrangian_engine.o: lagrangian_engine.cpp parameters.o velocity.o
	$(CC) -c lagrangian_engine.cpp -fopenmp

matrix_construction.o: matrix_construction.cpp parameters.o
	$(CC) -c matrix_construction.cpp

velocity.o: velocity.cpp velocity.h vectorXYZ.h date.h
	$(CC) -c velocity.cpp -lnetcdf_c++

ioutil.o: ioutil.cpp ioutil.h
	$(CC) -c ioutil.cpp 

VTK.o: VTK.cpp VTK.h
	$(CC) -c VTK.cpp

grid_construction.out: grid_construction.o parameters.o velocity.o VTK.o ioutil.o
	$(CC)  parameters.o grid_construction.o velocity.o VTK.o ioutil.o -o grid_construction.out $(LIBS)

lagrangian_engine.out: lagrangian_engine.o parameters.o velocity.o ioutil.o
	$(CC)  parameters.o lagrangian_engine.o velocity.o ioutil.o -o lagrangian_engine.out $(LIBS) -fopenmp

matrix_construction.out: matrix_construction.o parameters.o ioutil.o
	$(CC)  matrix_construction.o parameters.o ioutil.o -o matrix_construction.out

clean:
	$(RM) *.o 
