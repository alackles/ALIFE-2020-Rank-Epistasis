EMP_DIR := /home/austin/research/tools/empirical_clean/Empirical/source

CXX := g++-8

CFLAGS := -Wall -Wno-unused-function -iquote $(EMP_DIR)/ -std=c++17

OFLAGS_optim := -O3 -DNDEBUG
OFLAGS_debug := -g -pedantic -DEMP_TRACK_MEM  -Wnon-virtual-dtor -Wcast-align -Woverloaded-virtual

analysis: main.cc file_io.h organism.h nk.h edit_distance.h
	$(CXX) main.cc $(CFLAGS) $(OFLAGS_optim) -o analysis

debug: main.cc file_io.h organism.h nk.h edit_distance.h
	$(CXX) main.cc $(CFLAGS) $(OFLAGS_debug) -o analysis

clean:
	rm analysis
