ROOT:=../../
-include $(ROOT)/make/make_defs.mk

ALL_TARGETS:=prepare distgraph 

INCLUDES:=-I. -I../globals -I../sequential

all: $(ALL_TARGETS)

CXXFLAGS:=-fopenmp ${CXXFLAGS}

TARGET:=distgraph
SRCS:=graph_miner_mpi_omp_part.cpp main_mpi_omp_part.cpp ismin_mpi_omp_part.cpp external_neighbor_handler.cpp global_support_handler.cpp
LDFLAGS:=../globals/libglobals.a ../sequential/libseq2.a
ALL_DEPS:=$(ALL_DEPS) $(SRCS:.cpp=.dep)
-include $(ROOT)/make/compile_mpichexe.mk


-include $(ROOT)/make/make_tools.mk
-include $(addprefix deps/,$(ALL_DEPS))


