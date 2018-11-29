#!/bin/bash

# Runs QCE "single point" calculations.

# Get amf_ab as median of all optimized amf's (except for the pure ones).
get_amf_mix() {
    grep 'Best isobar' \
        {xw_0.187,xw_0.355,xw_0.522,xw_0.651,xw_0.742,xw_0.814,xw_0.944}/sampling.out | \
            awk '{
                a[NR]=$8
            }
            END {
                asort(a); 
                if (NR%2) {
                    print a[(NR+1)/2]
                } else {
                    print (a[NR/2]+a[NR/2+1])/2
                }
            }'
}

# Get amf_a and amf_b.
get_amf_a() {
    grep 'Best isobar' xw_1.0/sampling.out | awk '{print $7}'
}
get_amf_b() {
    grep 'Best isobar' xw_0.0/sampling.out | awk '{print $7}'
}

# Get bxv_a and bxv_b.
get_bxv_a() {
    grep 'Best isobar' xw_1.0/sampling.out | awk '{print $11}'
}
get_bxv_b() {
    grep 'Best isobar' xw_0.0/sampling.out | awk '{print $11}'
}

amf_ab=`get_amf_mix`
amf_a=`get_amf_a`
amf_b=`get_amf_b`
amf_pure="$amf_a $amf_b"
bxv_a=`get_bxv_a`
bxv_b=`get_bxv_b`
bxv_pure="$bxv_a $bxv_b"

# pure DMSO
echo xw_0.0
cd xw_0.0
    sed -i 's/    amf .*$/    amf '$amf_b'/' qce.input
    sed -i 's/    bxv .*$/    bxv '$bxv_b'/' qce.input
    ../../../peacemaker qce.input ../jcp.dmso.clusterset > qce.out
cd ..

# mixtures
for x in xw_0.187 xw_0.355 xw_0.522 xw_0.651 xw_0.742 xw_0.814 xw_0.944; do
    echo $x
    cd $x
        sed -i 's/    amf_pure .*$/    amf_pure '"$amf_pure"'/' qce.input
        sed -i 's/    amf .*$/    amf '$amf_ab'/' qce.input
        sed -i 's/    bxv_pure .*$/    bxv_pure '"$bxv_pure"'/' qce.input
        ../../../peacemaker qce.input ../jcp.mixed.clusterset > qce.out
    cd ..
done

# pure water
echo xw_1.0
cd xw_1.0
    sed -i 's/    amf .*$/    amf '$amf_a'/' qce.input
    sed -i 's/    bxv .*$/    bxv '$bxv_a'/' qce.input
    ../../../peacemaker qce.input ../jcp.water.clusterset > qce.out
cd ..
