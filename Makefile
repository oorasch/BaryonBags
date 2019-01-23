##################################################################################################
# makefile for a simple project
# available targets:
#      make
#      make clean
##################################################################################################
#  Edit the following lines according to your needs
##################################################################################################
PROGRAM	= main
SOURCES = main.cpp updates.cpp initialize.cpp debug.cpp parameters.cpp measure.cpp bags.cpp
OBJECTS = $(SOURCES:.cpp=.o)

CC	= gcc
CXX     = g++
F77	= gfortran

WARNINGS = -Wall -Weffc++ -Woverloaded-virtual -W -Wfloat-equal -Wshadow \
           -Wredundant-decls -Winline
#  -Wunreachable-code
# CXXFLAGS = -ffast-math -O3 -funroll-all-loops -DNDEBUG -msse3 ${WARNING} -ftree-vectorizer-verbose=2
#CXXFLAGS = -ffast-math -O3 -funroll-all-loops -DNDEBUG ${WARNING} -ftree-vectorizer-verbose=2 \
           # -fopenmp
CXXFLAGS = -std=c++11 -Ofast -march=native `gsl-config --cflags` -pedantic -Wall -Wextra -O3
#CXXFLAGS = -std=c++11
# CXXFLAGS = -g
#LDFLAGS  = -lcblas -fopenmp
LDFLAGS = `gsl-config --libs`
#LDFLAGS =  -lgslcblas -lm -lgsl

# CFLAGS	= -ffast-math -O3 -DNDEBUG -msse3 -fopenmp -fdump-tree-vect-details
# CFLAGS	= -ffast-math -O3 -funroll-loops -DNDEBUG -msse3 -fopenmp -ftree-vectorizer-verbose=2
# #CFLAGS	= -ffast-math -O3 -DNDEBUG -msse3 -fopenmp
# FFLAGS	= -ffast-math -O3 -DNDEBUG -msse3 -fopenmp
# LFLAGS  = -ffast-math -O3 -DNDEBUG -msse3 -fopenmp

##################################################################################################
#    Don't change anything below this line
##################################################################################################
all:	${PROGRAM}

${PROGRAM}:	${OBJECTS}
	$(CXX) $^ -o $@ ${LDFLAGS}

clean:
	rm -f ${PROGRAM} ${OBJECTS}

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<

.f.o:
	$(F77) -c $(FFLAGS) $<

###################################################################################################
#    some tools
# Cache behaviour (CXXFLAGS += -g  tracks down to source lines)
cache: ${PROGRAM}
	valgrind --tool=callgrind --simulate-cache=yes ./$^
#	kcachegrind callgrind.out.<pid> &

# Check for wrong memory accesses, memory leaks, ...
# use smaller data sets
mem: ${PROGRAM}
	valgrind -v --leak-check=yes --tool=memcheck --undef-value-errors=yes --track-origins=yes --log-file=$^.addr.out --show-reachable=yes ./$^

#  Simple run time profiling of your code
#  CXXFLAGS += -g -pg
#  LDFLAGS += -pg
prof: ${PROGRAM}
	./$^
	gprof -b ./$^ > gp.out
#	kprof -f gp.out -p gprof &
