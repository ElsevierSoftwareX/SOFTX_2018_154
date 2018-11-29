#!/bin/bash

# Calculates the Gibbs energies of mixing.

G0=`awk 'NR == 195 {print $3}' xw_0.0/thermo0.dat`
G1=`awk 'NR == 195 {print $3}' xw_1.0/thermo0.dat`


for x in 0.0 0.187 0.355 0.522 0.651 0.742 0.814 0.944 1.0; do
    G=`awk 'NR == 195 {print $3}' xw_$x/thermo0.dat`
    awk "BEGIN {print $x, $G-$x*$G1-(1.0-$x)*$G0}"
done
