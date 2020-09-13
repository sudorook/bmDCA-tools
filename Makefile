MSACOMPUTE := compute_msa_stats
ENERGY := compute_energies
DISTANCE := compare_distances
MSACOMPARE := compare_msa_stats
DMS := compute_dms
DCS := compute_dcs
HIST1D := compute_histogram
HIST2D := compute_histogram2d
AVERAGE := compute_average
SUBSET := subset_alignment
PATHFINDER := compute_paths
Z := compute_z

CC := g++

CXXFLAGS = -O3 $(shell pkg-config --cflags armadillo) \
           -fopenmp -std=c++11 -DARMA_DONT_USE_HDF5
LDFLAGS = -lm $(shell pkg-config --libs armadillo)

SRCPATH = src

SOURCES_MSACOMPUTE = ${SRCPATH}/compute_stats.cpp \
                     ${SRCPATH}/msa.cpp \
                     ${SRCPATH}/msa_stats.cpp \
                     ${SRCPATH}/utils.cpp
OBJECTS_MSACOMPUTE = $(SOURCES_MSACOMPUTE:%.cpp=%.o)

SOURCES_MSACOMPARE = ${SRCPATH}/compare_stats.cpp \
                     ${SRCPATH}/msa.cpp \
                     ${SRCPATH}/msa_compare.cpp \
                     ${SRCPATH}/utils.cpp
OBJECTS_MSACOMPARE = $(SOURCES_MSACOMPARE:%.cpp=%.o)

SOURCES_ENERGY = ${SRCPATH}/compute_energies.cpp \
                 ${SRCPATH}/msa.cpp \
                 ${SRCPATH}/utils.cpp
OBJECTS_ENERGY = $(SOURCES_ENERGY:%.cpp=%.o)

SOURCES_DISTANCE = ${SRCPATH}/compare_distances.cpp \
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

SOURCES_HIST1D = ${SRCPATH}/compute_histogram1d.cpp \
                 ${SRCPATH}/utils.cpp
OBJECTS_HIST1D = $(SOURCES_HIST1D:%.cpp=%.o)

SOURCES_HIST2D = ${SRCPATH}/compute_histogram2d.cpp \
                 ${SRCPATH}/utils.cpp
OBJECTS_HIST2D = $(SOURCES_HIST2D:%.cpp=%.o)

SOURCES_AVERAGE = ${SRCPATH}/average_models.cpp \
                  ${SRCPATH}/utils.cpp
OBJECTS_AVERAGE = $(SOURCES_AVERAGE:%.cpp=%.o)

SOURCES_SUBSET = ${SRCPATH}/subset_alignment.cpp \
                 ${SRCPATH}/msa.cpp \
                 ${SRCPATH}/utils.cpp
OBJECTS_SUBSET = $(SOURCES_SUBSET:%.cpp=%.o)

SOURCES_PATHFINDER = ${SRCPATH}/compute_paths.cpp \
                     ${SRCPATH}/msa.cpp \
                     ${SRCPATH}/utils.cpp
OBJECTS_PATHFINDER = $(SOURCES_PATHFINDER:%.cpp=%.o)

SOURCES_Z = ${SRCPATH}/compute_z.cpp \
            ${SRCPATH}/utils.cpp
OBJECTS_Z = $(SOURCES_Z:%.cpp=%.o)

.PHONY: all

all: $(MSACOMPUTE) $(MSACOMPARE) $(ENERGY) $(DISTANCE) $(DMS) $(DCS) $(HIST1D) $(HIST2D) $(AVERAGE) $(SUBSET) $(PATHFINDER)

$(MSACOMPUTE): $(OBJECTS_MSACOMPUTE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(MSACOMPARE): $(OBJECTS_MSACOMPARE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(ENERGY): $(OBJECTS_ENERGY)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(DISTANCE): $(OBJECTS_DISTANCE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(DMS): $(OBJECTS_DMS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(DCS): $(OBJECTS_DCS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(HIST1D): $(OBJECTS_HIST1D)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(HIST2D): $(OBJECTS_HIST2D)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(AVERAGE): $(OBJECTS_AVERAGE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(SUBSET): $(OBJECTS_SUBSET)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(PATHFINDER): $(OBJECTS_PATHFINDER)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(MSACOMPUTE)
	rm -f $(MSACOMPARE)
	rm -f $(ENERGY)
	rm -f $(DISTANCE)
	rm -f $(DMS)
	rm -f $(DCS)
	rm -f $(HIST1D)
	rm -f $(HIST2D)
	rm -f $(AVERAGE)
	rm -f $(SUBSET)
	rm -f $(PATHFINDER)
	rm -f $(OBJECTS_MSACOMPUTE)
	rm -f $(OBJECTS_MSACOMPARE)
	rm -f $(OBJECTS_ENERGY)
	rm -f $(OBJECTS_DISTANCE)
	rm -f $(OBJECTS_DMS)
	rm -f $(OBJECTS_DCS)
	rm -f $(OBJECTS_HIST1D)
	rm -f $(OBJECTS_HIST2D)
	rm -f $(OBJECTS_AVERAGE)
	rm -f $(OBJECTS_SUBSET)
	rm -f $(OBJECTS_PATHFINDER)
