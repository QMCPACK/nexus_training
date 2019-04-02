#! /usr/bin/env python

from nexus import obj
from nexus import read_structure
from nexus import generate_physical_system
from nexus import generate_pwscf_input

# generate using keywords + system

# read structure from file
s = read_structure('VO2_M1_afm.xsf')
s.elem[0] = 'V1' # set AFM pattern
s.elem[1] = 'V2'
s.elem[2] = 'V1'
s.elem[3] = 'V2'

# create physical system from structure
vo2 = generate_physical_system(
    structure = s,
    V1        = 13,
    V2        = 13,
    O         =  6,
    )

# generate pwscf input
pw = generate_pwscf_input(
    selector         = 'generic',
    calculation      = 'scf',
    disk_io          = 'low',
    verbosity        = 'high',
    wf_collect       = True,
    input_dft        = 'lda',
    hubbard_u        = obj(V1=3.5,V2=3.5),
    ecutwfc          = 350,
    bandfac          = 1.3,
    nosym            = True,
    occupations      = 'smearing',
    smearing         = 'fermi-dirac',
    degauss          = 0.0001,
    nspin            = 2,
    start_mag        = obj(V1=1.0,V2=-1.0),
    diagonalization  = 'david',
    conv_thr         = 1e-8,
    mixing_beta      = 0.2,
    electron_maxstep = 1000,
    system           = vo2,
    pseudos          = ['V.opt.upf','O.opt.upf'],
    kgrid            = (6,6,6),
    kshift           = (0,0,0),
    )

# manipulate and write
for ecut in [100,120,140,160]:
    pw.system.ecutwfc = ecut
    pw.write('scf_06_ecut_{0}.in'.format(ecut))
#end for
