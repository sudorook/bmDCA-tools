# SPDX-FileCopyrightText: 2020 - 2023 sudorook <daemon@nullcodon.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

MSACOMPUTE := compute_msa_stats
MSACOMPARE := compare_msa_stats
ENERGY := compute_energies
DISTANCE := compare_distances
GAUGE := compute_zero_gauge
DMS := compute_dms
DCS := compute_dcs
HIST1D := compute_histogram
HIST2D := compute_histogram2d
AVERAGE := compute_average
SUBSET := subset_alignment
PATHFINDER := compute_paths

CC := g++

CXXFLAGS = -O3 -Wall $(shell pkg-config --cflags armadillo) \
           -fopenmp -std=c++17 -DARMA_DONT_USE_HDF5 -DARMA_NO_DEBUG
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

SOURCES_GAUGE = ${SRCPATH}/compute_zero_gauge.cpp \
                ${SRCPATH}/utils.cpp
OBJECTS_GAUGE = $(SOURCES_GAUGE:%.cpp=%.o)

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

.PHONY: all

all: $(MSACOMPUTE) $(MSACOMPARE) $(ENERGY) $(DISTANCE) $(GAUGE) $(DMS) $(DCS) $(HIST1D) $(HIST2D) $(AVERAGE) $(SUBSET) $(PATHFINDER)

$(MSACOMPUTE): $(OBJECTS_MSACOMPUTE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(MSACOMPARE): $(OBJECTS_MSACOMPARE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(ENERGY): $(OBJECTS_ENERGY)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(DISTANCE): $(OBJECTS_DISTANCE)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(GAUGE): $(OBJECTS_GAUGE)
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
	rm -rf __pycache__
	rm -f $(MSACOMPUTE)
	rm -f $(MSACOMPARE)
	rm -f $(ENERGY)
	rm -f $(DISTANCE)
	rm -f $(GAUGE)
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
	rm -f $(OBJECTS_GAUGE)
	rm -f $(OBJECTS_DMS)
	rm -f $(OBJECTS_DCS)
	rm -f $(OBJECTS_HIST1D)
	rm -f $(OBJECTS_HIST2D)
	rm -f $(OBJECTS_AVERAGE)
	rm -f $(OBJECTS_SUBSET)
	rm -f $(OBJECTS_PATHFINDER)
