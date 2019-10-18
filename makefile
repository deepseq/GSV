CC = g++
IDIR = -Iinclude
LIBS = 

CPPFLAGS = $(IDIR) $(LIBS) -std=c++11 -c -Wall

SOURCES = src/main.cpp src/CmdExecution.cpp src/ProcessTool.cpp src/StringOperation.cpp src/CurrentTime.cpp src/DataPreparation.cpp src/MappingReads.cpp src/BaseCoverage.cpp src/MappingSummary.cpp src/ReadsAssembly.cpp src/BlastFilter.cpp src/CoverageFilter.cpp src/FilterSummary.cpp

OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = GSVMining

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm src/*.o $(EXECUTABLE)

