OBJS:=$(addprefix $(OBJSDIR)/,$(subst .cpp,.o,$(SRCS)))
DEPS:=$(addprefix $(DEPSDIR)/,$(subst .cpp,.dep,$(SRCS)))

#MPIOBJS:=$(addprefix $(OBJSDIR)/,$(subst .cpp,.o,$(MPI_SRCS)))
#MPIOBJS_$(TARGET):=$(MPIOBJS)
#MPI_TEST_OBJS_$(TARGET):=$(addprefix $(OBJSDIR)/,$(subst .cpp,.o,$(MPI_TEST_SRCS)))


OBJS_$(TARGET):=$(OBJS)
DEPS_$(TARGET):=$(TARGET_DEPS)
LDFLAGS_$(TARGET):=$(LDFLAGS)
CXXFLAGS_$(TARGET):=$(CXXFLAGS)
LD_DEPS_$(TARGET):=$(LD_DEPS)

$(TARGET):  $(OBJS_$(TARGET)) $(DEPS_$(TARGET)) $(LD_DEPS_$(TARGET))
	echo "LD $(@)"
	$(CXX)  -o $(@)  $(OBJS_$@) $(LDFLAGS_$@) 

