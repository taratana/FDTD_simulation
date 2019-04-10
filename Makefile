######################
#  Makefile Sample
#            project : name of an execution file generated finally
#            objs : series of program file name replaced ".cpp" with ".o"
#            CXX : compile command
#            CFLAGS : compile options
#
# How to use
#	make : compile and link
#	make clean : delete object files
#	make depend : analyze file interrelationship and generate a file (depend.inc) described dependence
######################

project = mato_FDTD_3D
objs= main.o FDTD.o GraphR.o
CXX=icpc
#CFLAGS= -O3 -ipo -xhost -mkl -Wall
CFLAGS= -O3 -ipo -qopenmp -xhost -mkl -Wall -std=c++11
#CFLAGS = -O3 -ipo -xS -openmp

.SUFFIXES: .cpp .o

$(project) : $(objs)
	$(CXX) -o $(project) $(CFLAGS)  $^

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(project) $(objs)

.PHONY: depend
depend: $(objs:.o=.cpp)
	-@ $(RM) depend.inc
	-@ for i in $^; do\
		cpp -MM $$i | sed "s/\ [_a-zA-Z0-9][_a-zA-Z0-9]*\.cpp//g" >> depend.inc;\
	done

-include depend.inc
