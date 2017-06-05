DEPSDIR:=deps
OBJSDIR:=objs

MPICXX:=mpicxx
CXX:=mpicxx
DEPCXX:=mpicxx
LD:=g++

CXXFLAGS:=-O4 -DLOG_INFO -DSTLPORT_UNAVAL -m64
DEPSCXXFLAGS:=-DLOG_TRACE
INCLUDES:=-I. -I$(ROOT)/src/globals
LDFLAGS:=


#GTEST_LIB:=$(ROOT)/gtest/distr_nostlp/lib/libgtest.a
#GTESTMAIN_LIB:=$(ROOT)/gtest/distr_nostlp/lib/libgtest_main.a
#GTEST_LD:=$(GTEST_LIB) $(GTESTMAIN_LIB)
#GTEST_INCLUDE:=-I$(ROOT)/gtest/distr_nostlp/include


