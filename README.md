# peakypaper

Code run to produce analyses in the manuscript "Fine mapping chromatin contacts in capture Hi-C data" by Eijsbouts et al. https://www.biorxiv.org/content/early/2018/01/05/243642

Note these commands make use of local software to run on our q (q.rb, qR.rb).  The details of these parts don't matter, you can extract the command that is being run.

## Sensitivity-specificity analysis in macrophage data

code in sens-spec

prepare input data (neg binom residuals) adjusted for bait2baits
```{sh}
cd sens-spec
for patt in N0 N1 N2 "'N0|N1'" "'N0|N2'" "'N1|N2'"; do
    qR.rb -r -c 6 prep-data.R --args patt="$patt"
done
```

generate list of baits for peaky - those with at least one small residual, FDR < 0.1, in a 10mb window around bait

```{sh}
qR.rb -r sens-spec/make-totest.R
```

## run each bait
```{sh}
./run.rb -r run
```

summarise output
```{sh}
qR.rb -a CWALLACE-SL2-CPU -p skylake-himem -t "04:00:00" collate-tests.R
cd ..
```

## Run peaky on CD4 T cell data

see cd4/Readme.sh

## Inference on peaky results
code to generate numbers, tables & figures for peaky paper all within the analysis directory.

```{sh}
./analysis/figure-NB-residuals-qqplots.R
./analysis/results-prepare-datafiles.R
./analysis/numbers-for-paper.R # fig 3
./analysis/assess.R ## table 2
./analysis/figures-examples.R ## figures 4,5
```

