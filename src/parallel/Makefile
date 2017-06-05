ROOT:=../../
-include $(ROOT)/make/make_defs.mk

ALL_TARGETS:=prepare pargraph_mpi_dyn pargraph_mpi_split_all pargraph_hybrid

INCLUDES:=-I. -I../sequential -I../globals

all: $(ALL_TARGETS)
CXXFLAGS:=-fopenmp ${CXXFLAGS}

TARGET:=pargraph_mpi_dyn
SRCS:=graph_miner_mpi_dyn.cpp main_mpi_dyn.cpp dynamic_load_balancing.cpp
LDFLAGS:=../globals/libglobals.a ../sequential/libseq.a
ALL_DEPS:=$(ALL_DEPS) $(SRCS:.cpp=.dep)
-include $(ROOT)/make/compile_mpichexe.mk

TARGET:=pargraph_mpi_split_all
SRCS:=graph_miner_mpi_split_all.cpp main_mpi_split_all.cpp dynamic_load_balancing.cpp
LDFLAGS:=../globals/libglobals.a ../sequential/libseq.a
ALL_DEPS:=$(ALL_DEPS) $(SRCS:.cpp=.dep)
-include $(ROOT)/make/compile_mpichexe.mk

TARGET:=pargraph_hybrid
SRCS:=graph_miner_mpi_omp_hybrid.cpp main_mpi_omp_hybrid.cpp ismin_mpi_omp_hybrid.cpp dynamic_hybrid_load_balancing_priv_params.cpp dynamic_hybrid_threads_load_balancing_priv_params.cpp
LDFLAGS:=../globals/libglobals.a ../sequential/libseq.a
ALL_DEPS:=$(ALL_DEPS) $(SRCS:.cpp=.dep)
-include $(ROOT)/make/compile_mpichexe.mk

-include $(ROOT)/make/make_tools.mk
-include $(addprefix deps/,$(ALL_DEPS))


