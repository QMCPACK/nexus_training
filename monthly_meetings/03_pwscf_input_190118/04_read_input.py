#! /usr/bin/env python

from pwscf_input import PwscfInput
from nexus import read_input

# read the input from a file
pw = PwscfInput('./Fe_start_ns_eig.in')


# reading this way also works
pw = read_input('./Fe_start_ns_eig.in',format='pwscf')


# manipulate and write
for ecut in [100,120,140,160]:
    pw.system.ecutwfc = ecut
    pw.write('scf_04_ecut_{0}.in'.format(ecut))
#end for
