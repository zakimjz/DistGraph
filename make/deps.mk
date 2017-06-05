deps/%.dep: %.cpp
	@echo DEP $<
	@mkdir -p deps
	@$(DEPCXX)   $(INCLUDES) $(DEPSCXXFLAGS) -M -MM -MT objs/$*.o  $^ > $@



deps/%.dep: tests/%.cpp
	@echo DEP $<
	@mkdir -p deps
	@$(DEPCXX)   $(INCLUDES) $(DEPSCXXFLAGS) -M  -MM -MT objs/$*.o $^ > $@


deps/%.dep: %.cu
	@echo DEP $<
	@mkdir -p deps
	@$(DEPCXX)   $(INCLUDES) $(CUDA_DEFS) -isystem /usr/local/cuda/include/  -M -MM -MT objs/$*.o   $^ > $@



-include $(shell ls deps/*.dep)
