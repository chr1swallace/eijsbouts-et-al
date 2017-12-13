# last file that should be created in each group
PREPDATA = dervied/actval-chicago.RData # last file created by results-prepare-datafiles.R
JOINDATA = derived/nonval-joined.RData
NUMBERS = derived/numbers.RData
FIGFM = figures/nonval-chic-mppc.pdf
FIGEX = figures/ex-IRF8-act-zoom1.pdf
FIGCOR = figures/MPPC-vs-outcome.pdf

# this is how we run Rscript
runscript = "Rscript $(1) > $(subst .R,.Rout,$(1))"

all : $(PREPDATA) $(JOINDATA) $(FIGFM) $(FIGEX) $(TABALL)
	echo "done"

## documented aliases
.PHONY: help prep joint finemap examples 

prep : $(PREPDATA) ## prepare data
	echo done
join : $(JOINDATA) ## join data to chicago and other annotations
	echo done
numbers : $(NUMBERS) ## final data for analysis, basic numbers and plots
	echo done
assess : $(FIGFM) ## assess by comparison to external measures
	echo done
exclear : ## remove existing example files
	rm figures/ex-*.pdf

## documented work
examples : figures-examples.R viz.R derived/nonprom-joined.RData ## make example figures
	q.rb -a tesla -r -c 1 -t '12:00:00' -j examples  $(call runscript,$<)

## the work
$(PREPDATA) : results-prepare-datafiles.R 
	q.rb -r -c 16 -t 02:00:00 -j prepdata $(call runscript,$<)

$(JOINDATA) : prepare-joint-data.R 
	q.rb -r -t 02:00:00 -j joindata $(call runscript,$<)

$(NUMBERS) : numbers-for-paper.R 
	q.rb -r -t 02:00:00 -j numbers $(call runscript,$<)

# $(NUMBERS) : numbers-for-paper.R
#         q.rb -r -t 01:00:00 -j numbers $(call runscript,$<)

$(FIGFM) : assess.R 
	q.rb -r -t  '01:00:00' -j finemap $(call runscript,$<)


# magic from https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
help:
	@grep -E '^[a-zA-Z_-]+ ?:.*?## .*$$' Makefile | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help

