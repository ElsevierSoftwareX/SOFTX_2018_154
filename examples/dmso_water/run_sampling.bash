#!/bin/bash

# Runs the QCE sampling calculations.

# Substitute amf and bxv ranges, if desired.
run_sed() {
    sed -i 's/    amf .*$/    amf 0.3 2.2 191/' sampling.input
    sed -i 's/    bxv .*$/    bxv 0.8 1.2 41/' sampling.input
}
run_sed2() {
    sed -i 's/    amf_mix .*$/    amf_mix 1.0 3.0 201/' sampling.input
    amf1=`awk '/Best isobar/ {print $7}' ../xw_1.0/sampling.out`
    amf2=`awk '/Best isobar/ {print $7}' ../xw_0.0/sampling.out`
    sed -i 's/    amf_pure .*$/    amf_pure '"$amf1 $amf2"'/' sampling.input
    bxv1=`awk '/Best isobar/ {print $11}' ../xw_1.0/sampling.out`
    bxv2=`awk '/Best isobar/ {print $11}' ../xw_0.0/sampling.out`
    sed -i 's/    bxv_pure .*$/    bxv_pure '"$bxv1 $bxv1"'/' sampling.input
}

# pure DMSO
echo xw_0.0
cd xw_0.0
    run_sed
    ../../../peacemaker sampling.input ../jcp.dmso.clusterset > sampling.out
cd ..

# pure water
echo xw_1.0
cd xw_1.0
    run_sed
    ../../../peacemaker sampling.input ../jcp.water.clusterset > sampling.out
cd ..

# mixtures
for x in xw_0.187 xw_0.355 xw_0.522 xw_0.651 xw_0.742 xw_0.814 xw_0.944; do
    echo $x
    cd $x
        run_sed2
        ../../../peacemaker sampling.input ../jcp.mixed.clusterset > sampling.out
    cd ..
done

