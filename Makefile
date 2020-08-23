COMPUTE := compute_msa_stats
ENERGY := compute_energies
DISTANCE := compute_distances
COMPARE := compare_msa_stats
DMS := compute_dms
DCS := compute_dcs
PATHFINDER := compute_paths

CC := g++

CXXFLAGS = -O3 $(shell pkg-config --cflags armadillo) \
           -fopenmp -std=c++11 -DARMA_DONT_USE_HDF5 -DARMA_NO_DEBUG
LDFLAGS = -lm $(shell pkg-config --libs armadillo)

SRCPATH = src

SOURCES_COMPUTE = ${SRCPATH}/compute_stats.cpp \
                  ${SRCPATH}/msa.cpp \
                  ${SRCPATH}/msa_stats.cpp \
                  ${SRCPATH}/utils.cpp
OBJECTS_COMPUTE = $(SOURCES_COMPUTE:%.cpp=%.o)

SOURCES_COMPARE = ${SRCPATH}/compare_stats.cpp \
                  ${SRCPATH}/msa.cpp \
                  ${SRCPATH}/msa_compare.cpp \
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

SOURCES_DMS = ${SRCPATH}/compute_dms.cpp \
              ${SRCPATH}/msa.cpp \
              ${SRCPATH}/msa_stats.cpp \
              ${SRCPATH}/utils.cpp
OBJECTS_DMS = $(SOURCES_DMS:%.cpp=%.o)

SOURCES_DCS = ${SRCPATH}/compute_dcs.cpp \
              ${SRCPATH}/msa.cpp \
              ${SRCPATH}/msa_stats.cpp \
              ${SRCPATH}/utils.cpp
OBJECTS_DCS = $(SOURCES_DCS:%.cpp=%.o)

SOURCES_PATHFINDER = ${SRCPATH}/compute_paths.cpp \
                     ${SRCPATH}/msa.cpp \
                     ${SRCPATH}/utils.cpp
OBJECTS_PATHFINDER = $(SOURCES_PATHFINDER:%.cpp=%.o)

.PHONY: all

all: $(COMPUTE) $(COMPARE) $(ENERGY) $(DISTANCE) $(DMS) $(DCS) $(PATHFINDER)

$(COMPUTE): $(OBJECTS_COMPUTE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(COMPARE): $(OBJECTS_COMPARE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(ENERGY): $(OBJECTS_ENERGY)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(DISTANCE): $(OBJECTS_DISTANCE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(DMS): $(OBJECTS_DMS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(DCS): $(OBJECTS_DCS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(PATHFINDER): $(OBJECTS_PATHFINDER)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(COMPUTE)
	rm -f $(COMPARE)
	rm -f $(ENERGY)
	rm -f $(DISTANCE)
	rm -f $(DMS)
	rm -f $(DCS)
	rm -f $(PATHFINDER)
	rm -f $(OBJECTS_COMPUTE)
	rm -f $(OBJECTS_COMPARE)
	rm -f $(OBJECTS_ENERGY)
	rm -f $(OBJECTS_DISTANCE)
	rm -f $(OBJECTS_DMS)
	rm -f $(OBJECTS_DCS)
	rm -f $(OBJECTS_PATHFINDER)
