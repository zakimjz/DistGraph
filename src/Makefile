#SILE:=-s
SILE:=-s
all: 
	@make $(SILE) -C globals $@ -j8
	@make $(SILE) -C sequential $@ -j8
	@make $(SILE) -C parallel $@ -j8
	@make $(SILE) -C distributed $@ -j8

deps: 
	@make $(SILE) -C globals $@ -j8
	@make $(SILE) -C sequential   $@ -j8
	@make $(SILE) -C parallel $@ -j8
	@make $(SILE) -C distributed $@ -j8


clean: 
	@make $(SILE) -C globals $@
	@make $(SILE) -C sequential $@
	@make $(SILE) -C parallel $@
	@make $(SILE) -C distributed $@
