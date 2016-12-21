#
# General C++ Makefile by Anton Crombach, a.b.m.crombach at uu.nl
#
# Program vars

MYPATH = /home/anton/local

##
## 64 bits
##

CXX = g++
# optimising flags (64bits)
CPPFLAGS = -Wall -O2 -march=x86-64 -mmmx -msse -msse2 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -D_LARGEFILE64_SOURCE 
#CPPFLAGS = -DDOUBLESTRANDBREAKS -Wall -O3 -march=x86-64 -mmmx -msse -msse2 -fprefetch-loop-arrays -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -D_LARGEFILE64_SOURCE

# debugging
#CPPFLAGS = -g -DDEBUG -Wall

# 'p' is for profiling information, ALSO needs pg in linker!
#CPPFLAGS = -Wall -O3 -march=x86-64 -pg
# linker options
#LNFLAGS = -pg

INCDIR = -I../include -I$(MYPATH)/include -I/usr/include
LIBDIR = -L$(MYPATH)/lib64 
LIBS = -lboost_program_options -lboost_filesystem -lboost_regex -lbgl-viz -lxerces-c 

# Source/object paths
vpath %.cc ../src ../test
vpath %.d .
vpath %.hh ../include
vpath %.o .

BINPATH = ../bin
LIBPATH = ../lib

# Project vars
PROJECT = blowhole64
LIBRARY = blowhole
ALL = distribution.o \
      main.o fluke.o config.o stream_manager.o model.o factory.o \
	  net_agent_reader.o net_pop_reader.o \
      observer_manager.o logger.o population.o environment.o \
      scaling.o selection.o \
      net_agent.o delta_agent.o simple_agent.o genreg_agent.o agent.o \
      genome.o chromosome.o repeat.o gene.o retroposon.o interaction.o \
	  hop_graph.o \
      observer.o subject.o
OBJECTS = $(ALL)


# Targets
all: $(PROJECT)

$(PROJECT): $(OBJECTS) 
	$(CXX) $(LNFLAGS) $(LIBDIR) $^ -o $(BINPATH)/$@ $(LIBS)

$(LIBRARY): $(OBJECTS)
	$(CXX) -shared -Wl,-soname,$@.so \
	$(LIBDIR) $(LIBS) $^ -o $(LIBPATH)/$@.so

$(OBJECTS): %.o: %.cc
	$(CXX) -c $(CPPFLAGS) $(INCDIR) $< -o $@

%.d: %.cc
	@set -e; rm -f $@; \
	$(CXX) $(INCDIR) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),realclean)
-include $(OBJECTS:.o=.d)
endif
endif

.PHONY: clean realclean 
clean:
	@rm -f *.d.* *.o *.d 

realclean: clean 
	@cd ../src; rm -f *~ .*.swp *.bak; \
	cd ../include; rm -f *~ .*.swp *.bak
