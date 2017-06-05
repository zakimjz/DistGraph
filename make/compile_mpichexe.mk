

OBJS:=$(addprefix $(OBJSDIR)/,$(subst .cpp,.o,$(SRCS)))
DEPS_$(TARGET):=$(addprefix $(DEPSDIR)/,$(subst .cpp,.dep,$(SRCS)))


CXX:=$(MPICXX)
EXETARGET_$(TARGET):=$(TARGET)

#EXETA_$(TARGET):=$(TARGET)
OBJS_$(TARGET):=$(addprefix $(OBJSDIR)/,$(subst .cpp,.o,$(SRCS)))
DEPS_$(TARGET):=$(DEPS)
LDFLAGS_$(TARGET):=$(LDFLAGS)
CXXFLAGS_$(TARGET):=$(CXXFLAGS)

MPIOBJS_$(TARGET):=$(addprefix $(OBJSDIR)/,$(subst .cpp,.o,$(MPI_SRCS)))
MPI_TEST_OBJS_$(TARGET):=$(addprefix $(OBJSDIR)/,$(subst .cpp,.o,$(MPI_TEST_SRCS)))
#MPIOBJS_$(TARGET):=$(MPIOBJS_$(TARGET)) $(MPI_TEST_OBJS_$(TARGET))


$(EXETARGET_$(TARGET)): $(OBJS_$(TARGET)) $(MPIOBJS_$(TARGET)) $(MPI_TEST_OBJS_$(TARGET))
	echo "LD MPICXX $(@)"
	$(MPICXX)  -o $@   $(CXXFLAGS_$(TARGET)) -Wl,--start-group  $(OBJS_$@) $(MPIOBJS_$@) $(MPI_TEST_OBJS_$@) $(LDFLAGS_$@) -Wl,--end-group


$(MPIOBJS_$(TARGET)):  $(MPI_SRCS)
	echo MPICXX $@
	$(MPICXX) -c $(subst .o,.cpp,$(@F)) -o $@ $(CXXFLAGS_$(TARGET))  $(INCLUDES)




$(MPI_TEST_OBJS_$(TARGET)): $(addprefix tests/,$(MPI_TEST_SRCS))
	echo MPICXX $@
	$(MPICXX) -c tests/$(subst .o,.cpp,$(@F)) -o $@ $(CXXFLAGS_$(TARGET)) $(INCLUDES)


