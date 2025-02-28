# Name of the final program
PROG := colloidgel

# Source files
SRC := main.cpp system.cpp utils.cpp mc.cpp

# Object files (same names, but with .o instead of .cpp)
OBJS := $(SRC:.cpp=.o)

# The C++ compiler
CXX := g++


# C++ compilation flags
CXXFLAGS := -std=c++11 -g  # add -O3, -Wall, etc. as needed
CXXFLAGS += -I/scratch/tacc/apps/gcc11_2/boost/1.86.0/include  # For Homebrew installations
CXXFLAGS += -I/scratch/tacc/apps/gcc11_2/gsl/2.8/include  # For Homebrew installations
# Linker flags (library search paths, RPATHs, etc.)
LDFLAGS := -L/scratch/tacc/apps/gcc11_2/gsl/2.8/lib
LDFLAGS += -L/scratch/tacc/apps/gcc11_2/boost/1.86.0/lib
# Libraries to link against
LDLIBS := -lboost_program_options -lgsl -lgslcblas -lm

# Default target
all: $(PROG)

# Link rule: build the final executable from object files
$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LDLIBS) -o $@

# Compile rule: build an object file from a .cpp
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Remove object files
clean:
	rm -f $(OBJS)

# Remove everything that can be regenerated + the final executable
distclean: clean
	rm -f $(PROG) *.debug *.txt *.log *.mol2 *.lammpstrj *.data

# Remove only runtime artifacts
rerunclean:
	rm -f *.debug *.txt *.log *.mol2 *.lammpstrj *.data
