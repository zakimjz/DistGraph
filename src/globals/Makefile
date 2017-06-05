ROOT:=../..
-include $(ROOT)/make/make_defs.mk

ALL_TARGETS:=libglobals.a #convert_to_partitions

all: prepare $(ALL_TARGETS)

TARGET:=libglobals.a
SRCS:=dbio.cpp \
      dfs_code.cpp \
      graph_output.cpp \
      graph_types.cpp \
      graph_getmin.cpp \
      logger.cpp \
      utils.cpp \
      memory_checker.cpp
ALL_DEPS:=$(ALL_DEPS) $(SRCS:.cpp=.dep)
-include $(ROOT)/make/compile_lib.mk

TARGET:=convert_to_partitions
SRCS:=convert_to_partitions.cpp
LDFLAGS:=libglobals.a
ALL_DEPS:=$(ALL_DEPS) $(SRCS:.cpp=.dep)
-include $(ROOT)/make/compile_mpichexe.mk

-include $(ROOT)/make/make_tools.mk

-include $(addprefix deps/,$(ALL_DEPS))

