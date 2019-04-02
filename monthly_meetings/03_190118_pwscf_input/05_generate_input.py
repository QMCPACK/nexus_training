#! /usr/bin/env python

from nexus import obj
from nexus import generate_pwscf_input

# generate using only keywords

pw = generate_pwscf_input(
    selector        = 'generic',
    # control inputs
    calculation     = 'scf',
    restart_mode    = 'from_scratch',
    wf_collect      = True,
    outdir          = './output',
    pseudo_dir      = '../pseudo/',
    prefix          = 'fe',
    etot_conv_thr   = 1.0e-9,
    forc_conv_thr   = 1.0e-6,
    tstress         = True,
    tprnfor         = True,
    # system inputs
    ibrav           = 1,
    celldm          = { 1 : 15 },
    nat             = 2,
    ntyp            = 1,
    ecutwfc         = 100,
    ecutrho         = 300,
    nbnd            = 18,
    occupations     = 'smearing',
    degauss         = 0.0005,
    smearing        = 'methfessel-paxton',
    nspin           = 2,
    assume_isolated = 'martyna-tuckerman',
    lda_plus_u      = True,
    hubbard_u       = { 1 : 3.1 },
    starting_magnetization = { 1 : 0.9 },
    starting_ns_eigenvalue = {(1,2,1) : 0.0,
                              (2,2,1) : 0.0476060,
                              (3,2,1) : 0.0476060,
                              (4,2,1) : 0.9654373,
                              (5,2,1) : 0.9954307},
    # electrons inputs
    conv_thr        = 1.0e-9,
    mixing_beta     = 0.7,
    diagonalization = 'david',
    mixing_fixed_ns = 500,
    # atomic_species inputs
    mass            = obj(Fe=58.69000),
    pseudos         = ['Fe.pbe-nd-rrkjus.UPF'],
    # atomic_positions inputs
    elem            = ['Fe','Fe'],
    pos             = [[2.070000000, 0.000000000, 0.000000000],    
                       [0.000000000, 0.000000000, 0.000000000]],
    pos_specifier   = 'angstrom',
    # k_points inputs
    kgrid           = (1,1,1),
    kshift          = (1,1,1),
    )


# manipulate and write
for ecut in [100,120,140,160]:
    pw.system.ecutwfc = ecut
    pw.write('scf_05_ecut_{0}.in'.format(ecut))
#end for
