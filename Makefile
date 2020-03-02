TARGET := compute_stats
CC := g++

CXXFLAGS = -O3 $(shell pkg-config --cflags armadillo)
LDFLAGS = -lm -fopenmp $(shell pkg-config --libs armadillo)

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
