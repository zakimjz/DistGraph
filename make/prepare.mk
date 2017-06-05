#.PHONY: prepare  objs deps

prepare: 
	mkdir -p objs
	mkdir -p deps

objs:
	mkdir -p objs

deps: $(wildcard *.dep)
	mkdir -p deps

clean:
	@rm -rf lib*.a *~ objs deps $(ALL_TARGETS)
	mkdir -p objs
	mkdir -p deps

.SILENT: all clean
