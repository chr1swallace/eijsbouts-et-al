## prepare input data
for patt in aCD4val aCD4pchic nCD4val nCD4pchic; do
    qR.rb -j prep$patt -r -t "6:00:00" -c 32 -p skylake-himem -a CWALLACE-SL2-CPU prep-data.R --args patt="$patt"
done

for patt in aCD4val aCD4pchic nCD4val nCD4pchic; do
    qR.rb -j mkt$patt -r -t "2:00:00" -c 32 -p skylake-himem -a CWALLACE-SL2-CPU make-totest.R --args patt="$patt"
done

## run things
./run.rb -r brun

## check correlation
./run.rb -r corr

## run reps 3/4
./run.rb -r brun

## additional fill ins
./run.rb -r run

## check correlation again
./run.rb -r corr

## make combined mppc.txt for analysis
./run.rb -r mppc

## count things
for f in output/*/rep-*; do n=`ls $f|wc -l`; echo $f $n; done
wc -w output/*/totest*

## copy output to derived
for f in output/*; do
    bf=`basename $f`
    cp $f/mppc.csv ~/bsu/peaky/derived/${bf}-mppc.csv
    cp $f/mppc-cor.txt ~/bsu/peaky/derived/${bf}-mppc-cor.txt
done

    
