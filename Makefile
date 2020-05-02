COMPUTE := compute_msa_stats
COMPARE := _compare_msa_stats
ENERGY := compute_energies
DISTANCE := compute_distances

CC := g++

CXXFLAGS = -O3 $(shell pkg-config --cflags armadillo) \
           -fopenmp -std=c++11 -DARMA_DONT_USE_HDF5
LDFLAGS = -lm $(shell pkg-config --libs armadillo)

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

SOURCES_DISTANCE = ${SRCPATH}/compute_distances.cpp \
                   ${SRCPATH}/msa.cpp \
                   ${SRCPATH}/utils.cpp
OBJECTS_DISTANCE = $(SOURCES_DISTANCE:%.cpp=%.o)

.PHONY: all

all: $(COMPUTE) $(COMPARE) $(ENERGY) $(DISTANCE)

$(COMPUTE): $(OBJECTS_COMPUTE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(COMPARE): $(OBJECTS_COMPARE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(ENERGY): $(OBJECTS_ENERGY)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(DISTANCE): $(OBJECTS_DISTANCE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(COMPUTE)
	rm -f $(COMPARE)
	rm -f $(ENERGY)
	rm -f $(DISTANCE)
	rm -f $(OBJECTS_COMPUTE)
	rm -f $(OBJECTS_COMPARE)
	rm -f $(OBJECTS_ENERGY)
	rm -f $(OBJECTS_DISTANCE)
