COMPUTE := compute_stats
COMPARE := compare_stats
ENERGY := compute_energies

CC := g++

CXXFLAGS = -O3 $(shell pkg-config --cflags armadillo)
LDFLAGS = -lm -fopenmp $(shell pkg-config --libs armadillo)

SRCPATH = src

SOURCES_COMPUTE = ${SRCPATH}/compute_stats.cpp \
                  ${SRCPATH}/msa.cpp \
                  ${SRCPATH}/msa_stats.cpp \
                  ${SRCPATH}/utils.cpp
OBJECTS_COMPUTE = $(SOURCES_COMPUTE:%.cpp=%.o)

SOURCES_COMPARE = ${SRCPATH}/compare_stats.cpp \
                  ${SRCPATH}/msa.cpp \
                  ${SRCPATH}/msa_stats.cpp \
                  ${SRCPATH}/utils.cpp
OBJECTS_COMPARE = $(SOURCES_COMPARE:%.cpp=%.o)

SOURCES_ENERGY = ${SRCPATH}/compute_energies.cpp \
                 ${SRCPATH}/msa.cpp \
                 ${SRCPATH}/utils.cpp
OBJECTS_ENERGY = $(SOURCES_ENERGY:%.cpp=%.o)

.PHONY: all

all: $(COMPUTE) $(COMPARE) $(ENERGY)

$(COMPUTE): $(OBJECTS_COMPUTE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(COMPARE): $(OBJECTS_COMPARE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(ENERGY): $(OBJECTS_ENERGY)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(COMPUTE)
	rm -f $(COMPARE)
	rm -f $(ENERGY)
	rm -f $(OBJECTS_COMPUTE)
	rm -f $(OBJECTS_COMPARE)
	rm -f $(OBJECTS_ENERGY)
