JULIA = julia
MAIN = range_funcs.jl
INPUT ?= input.txt
METHOD1 ?= Lagrange3
METHOD2 ?= Taylor4

.PHONY: intervals analyze compare default

default: intervals

intervals:
	$(JULIA) $(MAIN) intervals $(INPUT)

analyze:
	$(JULIA) $(MAIN) analyze $(METHOD1) $(METHOD2) $(INPUT)

compare:
	$(JULIA) $(MAIN) compare $(METHOD1) $(METHOD2) $(INPUT)
