# Name of the final program
PROG := colloidgel

# Directories
SRC_DIR := src
INCLUDE_DIR := include
BUILD_DIR := build
BIN_DIR := bin
# Installation directories
PREFIX := /usr/local
BINDIR := $(PREFIX)/bin
INCLUDEDIR := $(PREFIX)/include/$(PROG)
# Source files (all .cpp files in SRC_DIR)
SRC := $(wildcard $(SRC_DIR)/*.cpp)

# Object files (same names, but with .o instead of .cpp and in BUILD_DIR)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC))

# The C++ compiler
CXX := clang++

# C++ compilation flags
CXXFLAGS := -std=c++11 -g -Wall -O3
CXXFLAGS += -I$(INCLUDE_DIR)
CXXFLAGS += -I/opt/homebrew/include  # For Homebrew installations

# Linker flags (library search paths, RPATHs, etc.)
LDFLAGS := -L/opt/homebrew/lib
LDFLAGS += -L/usr/local/lib

# Libraries to link against
LDLIBS := -lboost_program_options -lgsl -lgslcblas -lm

# Default target
all: directories $(BIN_DIR)/$(PROG)

# Create necessary directories
directories:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)

# Link rule: build the final executable from object files
$(BIN_DIR)/$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LDLIBS) -o $@
	@echo "Build complete: $(BIN_DIR)/$(PROG)"

# Compile rule: build an object file from a .cpp
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Remove object files
clean:
	rm -rf $(BUILD_DIR)

# Remove everything that can be regenerated + the final executable
distclean: clean
	rm -rf $(BIN_DIR)
	rm -f *.debug *.txt *.log *.mol2 *.lammpstrj *.data

# Remove only runtime artifacts
rerunclean:
	rm -f *.debug *.txt *.log *.mol2 *.lammpstrj *.data

# Phony targets
.PHONY: all clean distclean rerunclean directories

# Dependencies
-include $(OBJS:.o=.d)

# Generate dependency files
$(BUILD_DIR)/%.d: $(SRC_DIR)/%.cpp
	@set -e; rm -f $@; \
	$(CXX) -MM -MP -MT '$(@:.d=.o)' $(CXXFLAGS) $< > $@

# Install the executable and header files to system directories
install: $(BIN_DIR)/$(PROG)
	@echo "Installing $(BIN_DIR)/$(PROG) to $(BINDIR)..."
	@sudo mkdir -p $(BINDIR)
	@sudo install -m 755 $(BIN_DIR)/$(PROG) $(BINDIR)
	@echo "Installing header files to $(INCLUDEDIR)..."
	@sudo mkdir -p $(INCLUDEDIR)
	@sudo install -m 644 include/*.h $(INCLUDEDIR)
	@echo "Installation complete."

# Uninstall the program and header files
uninstall:
	@echo "Removing $(PROG) from $(BINDIR)..."
	@sudo rm -f $(BINDIR)/$(PROG)
	@echo "Removing header files from $(INCLUDEDIR)..."
	@sudo rm -rf $(INCLUDEDIR)
	@echo "Uninstallation complete."

# Help target
help:
	@echo "Available targets:"
	@echo "  all        - Build the program (default)"
	@echo "  install    - Install the executable to $(BINDIR) and headers to $(INCLUDEDIR)"
	@echo "  uninstall  - Remove the installed executable and header files"
	@echo "  clean      - Remove object files"
	@echo "  distclean  - Remove object files, executable, and all generated files"
	@echo "  rerunclean - Remove only runtime artifacts"
	@echo "  help       - Display this help message"