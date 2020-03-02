TARGET := compute_stats
CC := g++

CXXFLAGS = -Wall -Wextra -O3
LDFLAGS = -lm -fopenmp -larmadillo

SRCPATH = src

# HEADERS = $(shell find $(SRCPATH) -name '*.hpp' | sort -k 1nr | cut -f2-)
SOURCES = $(shell find $(SRCPATH) -name '*.cpp' | sort -k 1nr | cut -f2-)
OBJECTS = $(SOURCES:%.cpp=%.o)

.PHONY: all

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(TARGET)
	rm -f $(OBJECTS)
