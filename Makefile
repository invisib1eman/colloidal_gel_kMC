PROG1 = COLLOIDGEL

SRC1 = main.cpp system.cpp utils.cpp mc.cpp

OBJS1 = ${SRC1:.cpp=.o}

CXX = clang++
# Compilation flags
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -I/opt/homebrew/include

# Linker flags
LDFLAGS = -L/opt/homebrew/lib -lboost_program_options -lgsl -lgslcblas -lm


 
#CXXFLAGS=-O3 -funroll-loops -DNDEBUG 


all: $(PROG1) 

$(PROG1):  $(OBJS1)
	 $(CXX) -o $@ $^ $(LDFLAGS)

%.o:  %.cpp
	$(CXX) -c $< $(CXXFLAGS)

clean: 
	rm -rf *.o

distclean:
	rm -f $(PROG1) *.o *.debug *.txt *.log *.mol2 *.lammpstrj

rerunclean:
	rm -f *.debug *.txt *.log *.mol2 *.lammpstrj
