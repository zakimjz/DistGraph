objs/%.o: tests/%.cpp  deps/%.dep
	@echo $(CXX) $<
	@$(CXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDES)



objs/%.o: %.cpp deps/%.dep
	@echo $(CXX) $<
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@



objs/%.o: %.cu deps/%.dep
	@echo NVCC $<
	@$(CXX) $(CXXFLAGS)  $(INCLUDES) -c $< -o $@ 



