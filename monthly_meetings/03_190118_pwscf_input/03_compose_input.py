#! /usr/bin/env python

from numpy import array
from generic import obj
from pwscf_input import PwscfInput

# directly compose input in Nexus internal representation

pw = PwscfInput()
pw.control.set(
    calculation   = 'scf',
    restart_mode  = 'from_scratch',
    wf_collect    = True,
    outdir        = './output',
    pseudo_dir    = '../pseudo/',
    prefix        = 'fe',
    etot_conv_thr = 1.0e-9,
    forc_conv_thr = 1.0e-6,
    tstress       = True,
    tprnfor       = True,
    )
pw.system.set(
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
    starting_ns_eigenvalue = { (1,2,1) : 0.0,
                               (2,2,1) : 0.0476060,
                               (3,2,1) : 0.0476060,
                               (4,2,1) : 0.9654373,
                               (5,2,1) : 0.9954307},
    )
pw.electrons.set(
    conv_thr        = 1.0e-9,
    mixing_beta     = 0.7,
    diagonalization = 'david',
    mixing_fixed_ns = 500,
    )
pw.atomic_species.set(
    atoms            = ['Fe'],
    masses           = obj(Fe=58.69000),
    pseudopotentials = obj(Fe='Fe.pbe-nd-rrkjus.UPF'),
    )
pw.atomic_positions.set(
    specifier = 'angstrom',
    atoms     = ['Fe','Fe'],
    positions = array([
        [2.070000000,   0.000000000,   0.000000000],   
        [0.000000000,   0.000000000,   0.000000000], 
        ]),
    )
pw.k_points.set(
    specifier = 'automatic',
    grid      = (1,1,1),
    shift     = (1,1,1),
    )


# print internal representation
print 
print 'internal representation'
print pw  # this works for all PwscfInput's, however obtained

# manipulate and write
for ecut in [100,120,140,160]:
    pw.system.ecutwfc = ecut
    pw.write('scf_03_ecut_{0}.in'.format(ecut))
#end for



# this way also works
#
#pw = PwscfInput()
#pw.control.set(
#    calculation   = 'scf' ,
#    restart_mode  = 'from_scratch' ,
#    wf_collect    = True ,
#    outdir        = './output' ,
#    pseudo_dir    = '../pseudo/' ,
#    prefix        = 'fe' ,
#    etot_conv_thr = 1.0e-9 ,
#    forc_conv_thr = 1.0e-6 ,
#    tstress       = True ,
#    tprnfor       = True ,
#    )
#pw.system.set(
#    ibrav           = 1,
#    nat             = 2,
#    ntyp            = 1,
#    ecutwfc         = 100 ,
#    ecutrho         = 300 ,
#    nbnd            = 18,
#    occupations     = 'smearing',
#    degauss         = 0.0005 ,
#    smearing        = 'methfessel-paxton' ,
#    nspin           = 2 ,
#    assume_isolated = 'martyna-tuckerman',
#    lda_plus_u      = True ,
#    )
#pw.system.set({
#    'celldm(1)' : 15,
#    'starting_magnetization(1)' : 0.9,
#    'hubbard_u(1)' : 3.1,
#    'starting_ns_eigenvalue(1,2,1)' : 0.0,
#    'starting_ns_eigenvalue(2,2,1)' : 0.0476060,
#    'starting_ns_eigenvalue(3,2,1)' : 0.0476060,
#    'starting_ns_eigenvalue(4,2,1)' : 0.9654373,
#    'starting_ns_eigenvalue(5,2,1)' : 0.9954307,
#    })
#pw.electrons.set(
#    conv_thr        = 1.0e-9 ,
#    mixing_beta     = 0.7 ,
#    diagonalization = 'david' ,
#    mixing_fixed_ns = 500,
#    )
#pw.atomic_species.set(
#    atoms            = ['Fe'],
#    masses           = obj(Fe=58.69000),
#    pseudopotentials = obj(Fe='Fe.pbe-nd-rrkjus.UPF'),
#    )
#pw.atomic_positions.set(
#    specifier = 'angstrom',
#    atoms     = ['Fe','Fe'],
#    positions = array([
#        [2.070000000,   0.000000000,   0.000000000],   
#        [0.000000000,   0.000000000,   0.000000000], 
#        ]),
#    )
#pw.k_points.set(
#    specifier = 'automatic',
#    grid      = (1,1,1),
#    shift     = (1,1,1),
#    )
