#!/usr/bin/env python

#=========================================================================================
# Peacemaker -- A Quantum Cluster Equilibrium Code.
#
# Copyright 2004-2006 Barbara Kirchner, University of Bonn
# Copyright 2007-2012 Barbara Kirchner, University of Leipzig
# Copyright 2013-2016 Barbara Kirchner, University of Bonn
#
# This file is part of Peacemaker.
#
# Peacemaker is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Peacemaker is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Peacemaker.  If not, see <http://www.gnu.org/licenses/>
#=========================================================================================

import sys
import os.path
import math

# Values taken from Bondi's compilation.

radii = {"H"  : 1.20,
         "He" : 1.40,
         "C"  : 1.70,
         "N"  : 1.55,
         "O"  : 1.52,
         "F"  : 1.47,
         "Ne" : 1.54,
         "Si" : 2.10,
         "P"  : 1.80,
         "S"  : 1.80,
         "Cl" : 1.75,
         "Ar" : 1.88,
         "As" : 1.85,
         "Se" : 1.90,
         "Br" : 1.85,
         "Kr" : 2.02,
         "Te" : 2.06,
         "I"  : 1.98,
         "Xe" : 2.16}

if len(sys.argv) != 2:
    print("Usage: {} <file>".format(sys.argv[0]))
    sys.exit(1)

filename = sys.argv[1]

with open(filename, "r") as fin:
    line = fin.readline()
    print(line)
    natoms = int(line)
    line = fin.readline()
    volume = 0.0
    for i in range(natoms):
        line = fin.readline()
        label, x, y, z = line.split()
        r = radii[label]
        print(x, y, z, r)
