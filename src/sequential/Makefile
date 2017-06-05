ROOT:=../../
-include $(ROOT)/make/make_defs.mk

ALL_TARGETS:=prepare graph_miner_seq libseq.a graph_miner_seq2 libseq2.a

INCLUDES:=-I. -I../globals

all: $(ALL_TARGETS)


TARGET:=graph_miner_seq
SRCS:=graph_miner.cpp main.cpp ismin.cpp
LDFLAGS:=../globals/libglobals.a 
ALL_DEPS:=$(ALL_DEPS) $(SRCS:.cpp=.dep)
-include $(ROOT)/make/compile_mpichexe.mk

TARGET:=libseq.a
SRCS:=graph_miner.cpp ismin.cpp
LDFLAGS:=../globals/libglobals.a 
ALL_DEPS:=$(ALL_DEPS) $(SRCS:.cpp=.dep)
-include $(ROOT)/make/compile_lib.mk

TARGET:=graph_miner_seq2
SRCS:=graph_miner2.cpp main2.cpp ismin2.cpp
LDFLAGS:=../globals/libglobals.a 
ALL_DEPS:=$(ALL_DEPS) $(SRCS:.cpp=.dep)
-include $(ROOT)/make/compile_mpichexe.mk

TARGET:=libseq2.a
SRCS:=graph_miner2.cpp ismin2.cpp
LDFLAGS:=../globals/libglobals.a 
ALL_DEPS:=$(ALL_DEPS) $(SRCS:.cpp=.dep)
-include $(ROOT)/make/compile_lib.mk

-include $(ROOT)/make/make_tools.mk
-include $(addprefix deps/,$(ALL_DEPS))


