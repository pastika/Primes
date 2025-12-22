# Define the C++ compiler
CXX = g++

# Define compiler flags (optional)
CXXFLAGS = -O2 -g

# Find all .cpp files in the current directory
# and use patsubst to generate the desired executable names
SRCS := $(wildcard *.cpp)
EXES := $(basename $(SRCS))

.PHONY: all clean

# Default target: builds all executables
all: $(EXES)

# Phony target to clean up generated files
clean:
	rm -f $(EXES) *.o

# Generic rule to compile any executable from its corresponding .cpp file
# $@ is the target name (executable name)
# $< is the first prerequisite (source file)
%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@
