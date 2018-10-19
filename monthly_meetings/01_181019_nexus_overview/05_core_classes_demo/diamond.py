#! /usr/bin/env python

from nexus import settings,job,run_project
from nexus import read_structure,get_machine
from nexus import generate_physical_system
from nexus import generate_pwscf

b = '\n'+40*'='+'\n'

settings(
    pseudo_dir  = './pseudopotentials',
    machine     = 'ws16'
    #machine     = 'titan',
    #account     = 'XYZ123',
    )

if settings.machine=='ws16':
    scf_job = job(cores=16,app='pw.x')
elif settings.machine=='titan':
    titan = get_machine('titan')
    titan.queue_size = 2
    scf_job = job(nodes=1,hours=1,app='pw.x')
else:
    print 'no jobs for machine!'
#end if

s = read_structure('d16bulk.POSCAR')
s.remove([[0.8925,0.8925,0.8925]]) 
print b+'structure w/ vacancy\n',s
s.write('d16vac.xsf')

dia16vac = generate_physical_system(
    structure = s,
    C         = 4,
    )
print b+'net charge:',dia16vac.net_charge
              
scf = generate_pwscf(
    identifier  = 'scf',
    path        = 'diamond/scf',
    job         = scf_job,
    input_type  = 'generic',
    calculation = 'scf',
    input_dft   = 'lda', 
    ecutwfc     = 200,
    mixing_beta = 0.2,
    degauss     = 0.01,
    conv_thr    = 1e-6,
    system      = dia16vac,
    kgrid       = (3,3,3),
    kshift      = (0,0,0),
    pseudos     = ['C.BFD.upf'], 
    )
print b+'scf directory:',scf.locdir

si = scf.input
print b+'scf input data\n',si
print b+'scf input file\n',si.write()

pm = run_project()
print b+'status after execute:'
pm.write_simulation_status()

sa = scf.load_analyzer_image()
print b+'DFT Energy:',sa.E
