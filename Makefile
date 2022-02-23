# Various flags
#CXX  = clang++
CXX = g++
LINK = $(CXX)
#CXXFLAGS = -std=c++11 -Wall -g 
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -pedantic -Wno-unknown-pragmas
#CXXFLAGS = -std=c++11 -Wall -O3 -fopenmp
LFLAGS = -lm
#LFLAGS = -lm -fopenmp
# ifneq "$(findstring noomp, $(MAKECMDGOALS))" "noomp"
# 	CXXFLAGS += -fopenmp
# 	LFLAGS += -fopenmp
# endif

TARGET = partition-validation

HEADER = partition-validation.h
FILES = partition-validation.cc

OBJECTS = $(FILES:.cc=.o)

$(TARGET): ${OBJECTS}
	$(LINK) $^ $(LFLAGS) -o $@

all: $(TARGET)

clean:
	rm -f $(OBJECTS)

distclean:
	rm -f $(OBJECTS) $(TARGET)

# Compile and dependency
$(OBJECTS): $(HEADER) Makefile



