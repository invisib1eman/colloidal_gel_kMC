<<<<<<< HEAD
# Name of the final program
PROG := NANOROD

# Source files
SRC := main.cpp system.cpp utils.cpp mc.cpp quarternion.cpp
=======
PROG1 = COLLOIDGEL

SRC1 = main.cpp system.cpp utils.cpp mc.cpp
>>>>>>> origin/main

# Object files (same names, but with .o instead of .cpp)
OBJS := $(SRC:.cpp=.o)

<<<<<<< HEAD
# The C++ compiler
CXX := clang++
=======
CXX = clang++
# Compilation flags
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -I/opt/homebrew/include

# Linker flags
LDFLAGS = -L/opt/homebrew/lib -lboost_program_options -lgsl -lgslcblas -lm
>>>>>>> origin/main

# Preprocessor flags (includes)
CPPFLAGS := -I/opt/homebrew/include

# C++ compilation flags
CXXFLAGS := -std=c++11 -g  # add -O3, -Wall, etc. as needed
CXXFLAGS += -I/opt/homebrew/include  # For Homebrew installations

# Linker flags (library search paths, RPATHs, etc.)
LDFLAGS := -L/opt/homebrew/lib
LDFLAGS += -L/usr/local/lib -lgsl -lgslcblas -lboost_program_options

# Libraries to link against
LDLIBS := -lboost_program_options -lgsl -lgslcblas -lm

<<<<<<< HEAD
# Default target
all: $(PROG)

# Link rule: build the final executable from object files
$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LDLIBS) -o $@
=======
$(PROG1):  $(OBJS1)
	 $(CXX) -o $@ $^ $(LDFLAGS)

%.o:  %.cpp
	$(CXX) -c $< $(CXXFLAGS)
>>>>>>> origin/main

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