PROG1 = NANOROD

SRC1 = main.cpp system.cpp utils.cpp mc.cpp quarternion.cpp

OBJS1 = ${SRC1:.cpp=.o}

CXX = g++
CXXFLAGS= -lboost_program_options -lgsl -lgslcblas -std=c++11 -lm
#LIBS=-I/usr/local/gsl/include/
#on Mac M1
LIBS=-I/opt/homebrew/include -L/opt/homebrew/lib
#-L/usr/local/gsl/lib


 
#CXXFLAGS=-O3 -funroll-loops -DNDEBUG 


all: $(PROG1) 

$(PROG1):  $(OBJS1)
	 $(CXX) $(LIBS)  $(CXXFLAGS) -o $@ $^

%.o:  %.cpp
	$(CXX) $(LIBS) -c -o $@ $<

clean: 
	rm -rf *.o 

distclean:
	rm -f $(PROG1) *.o *.debug
