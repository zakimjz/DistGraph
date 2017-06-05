# OBJS, COMPILER, TARGET, FLAGS, DEPS, CXX, AR, INCLUDES, CXXFLAGS

OBJS:=$(addprefix $(OBJSDIR)/,$(subst .cpp,.o,$(SRCS)))
DEPS:=$(addprefix $(DEPSDIR)/,$(subst .cpp,.dep,$(SRCS)))

#MPI_SRCS_$(TARGET):=$(MPI_SRCS)

#MPIOBJS:=$(addprefix $(OBJSDIR)/,$(subst .cpp,.o,$(MPI_SRCS)))
#MPIOBJS_$(TARGET):=$(MPIOBJS)

TARGET_$(TARGET):=$(TARGET)
OBJS_$(TARGET):=$(addprefix $(OBJSDIR)/,$(subst .cpp,.o,$(SRCS)))
CXXFLAGS_$(TARGET):=$(CXXFLAGS)

$(TARGET): $(OBJS_$(TARGET)) 
	@echo AR $@
	$(AR) rcs $@ $(OBJS_$@)


#$(MPIOBJS_$(TARGET)): $(MPI_SRCS_$(TARGET))
#	echo MPICXX $@ 
#	$(MPICXX) -c $(subst .o,.cpp,$(@F)) -o $@ $(CXXFLAGS_$(TARGET)) $(INCLUDES)



