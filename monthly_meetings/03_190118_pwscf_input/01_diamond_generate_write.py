#! /usr/bin/env python

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pwscf_input

settings(
    pseudo_dir    = './pseudopotentials',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws16'
    )

dia16 = generate_physical_system(
    units  = 'A',
    axes   = [[ 1.785,  1.785,  0.   ],
              [ 0.   ,  1.785,  1.785],
              [ 1.785,  0.   ,  1.785]],
    elem   = ['C','C'],
    pos    = [[ 0.    ,  0.    ,  0.    ],
              [ 0.8925,  0.8925,  0.8925]],
    tiling = (2,2,2),
    kgrid  = (1,1,1),
    kshift = (0,0,0),
    C      = 4
    )

# most familiar
scf = generate_pwscf(
    # Nexus Simulation class inputs
    identifier   = 'scf',
    path         = 'diamond/scf',
    job          = job(cores=16,app='pw.x'),
    # PwscfInput class inputs
    input_type   = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8, 
    nosym        = True,
    wf_collect   = True,
    system       = dia16,         # shared by Simulation/PwscfInput
    pseudos      = ['C.BFD.upf'], # shared by Simulation/PwscfInput
    )

print
print '===================='
print 'overview of contents'
print '===================='
print repr(scf.input) # PwscfInput class


#equivalent to the above
scf_input = generate_pwscf_input(
    # PwscfInput class inputs
    selector     = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8, 
    nosym        = True,
    wf_collect   = True,
    system       = dia16,
    pseudos      = ['C.BFD.upf'], 
    )

scf = generate_pwscf(
    # Nexus Simulation class inputs
    identifier   = 'scf',
    path         = 'diamond/scf2',
    job          = job(cores=16,app='pw.x'),
    pseudos      = ['C.BFD.upf'], # pseudos also here for file transfer
    system       = dia16,         # system info cached in sim object
    input        = scf_input,
    )

print
print '============================='
print 'actually write the input file'
print '============================='
text = scf_input.write() # write to string
print text
scf_input.write('./scf_01.in') # write to file



# no need to run, input demo
#run_project()
