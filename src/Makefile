ifneq ($(MAKECMDGOALS),clean)
  # test if R was built as a library
  RCMDERR := $(shell R CMD config --cppflags 2>&1 > /dev/null)

  ifeq ($(RCMDERR),)
    CPPFLAGS := $(shell R CMD config --cppflags)
    LDFLAGS := $(shell R CMD config --ldflags)
  else
    CPPFLAGS = -DMATHLIB_STANDALONE
    LDFLAGS = -lRmath
    ifndef MAKE_RESTARTS
      $(warning $(RCMDERR))
      $(info )
      $(info R was not built as a library, or was not found at all.)
      $(info Trying to use the standalone R math library (libRmath) instead.)
      $(info If the compilation fails, refer to the online documentation.)
      $(info )
    endif
  endif
endif

CXXFLAGS = -O3 -Wall -pedantic
LDLIBS = -lz
PROGRAMS = taco
SOURCES = $(wildcard *.cpp)

all: $(PROGRAMS)

taco : taco.o Genome.o HypothesesSet.o NarrowPeak.o PositionWeightMatrix.o Specification.o StatsSet.o
	$(CXX) -o $@ $^ -pthread $(LDFLAGS) $(LDLIBS)

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CXX) $(CPPFLAGS) -MM -o $@.$$$$ $<; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

ifneq ($(MAKECMDGOALS),clean)
  -include $(SOURCES:.cpp=.d)
endif

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

clean:
	rm -f $(PROGRAMS) *.o *.d

.PHONY: all clean
