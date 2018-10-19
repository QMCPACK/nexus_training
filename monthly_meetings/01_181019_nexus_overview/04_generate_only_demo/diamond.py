#! /usr/bin/env python

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack,vmc

settings(
    pseudo_dir    = './pseudopotentials',
    status_only   = 0,
    generate_only = 1,
    sleep         = .3,
    machine       = 'oic5'
    )

if settings.machine=='ws16':
    scf_job  = job(cores=16,app='pw.x')
    conv_job = job(cores=1,app='pw2qmcpack.x')
    qmc_job  = job(cores=16,threads=4,app='qmcpack')
elif settings.machine=='oic5':
    pwscf_modules = '''
module ()
{
    eval `/opt/modules/3.1.6/bin/modulecmd bash $*`
}
module purge
module load mpi/openmpi-1.4.5-pgi
module load PGI/2013-64bit
module load composerxe/2013.5.192
module load hdf5/1.8.8-pgi-parallel
'''
    qmcpack_modules = '''
module ()
{
    eval `/opt/modules/3.1.6/bin/modulecmd bash $*`
}
module purge
module load composerxe/2013.5.192
module load mpi/openmpi-1.6.4-gcc4
'''
    scf_job  = job(nodes=1,minutes=30,presub=pwscf_modules,app='pw.x')
    conv_job = job(cores=1,minutes=30,presub=pwscf_modules,app='pw2qmcpack.x')
    qmc_job  = job(nodes=1,threads=4,minutes=30,presub=qmcpack_modules,app='qmcpack')
else:
    print 'no jobs for machine!'
    exit()
#end if


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
              
scf = generate_pwscf(
    identifier   = 'scf',
    path         = 'diamond/scf',
    job          = scf_job,
    input_type   = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8, 
    nosym        = True,
    wf_collect   = True,
    system       = dia16,
    pseudos      = ['C.BFD.upf'], 
    )

conv = generate_pw2qmcpack(
    identifier   = 'conv',
    path         = 'diamond/scf',
    job          = conv_job,
    write_psir   = False,
    dependencies = (scf,'other'),
    )

qmc = generate_qmcpack(
    identifier   = 'vmc',
    path         = 'diamond/vmc',
    job          = qmc_job,
    input_type   = 'basic',
    system       = dia16,
    pseudos      = ['C.BFD.xml'],
    orbitals_h5  = '../scf/pwscf_output/pwscf.pwscf.h5',
    jastrows     = [],
    calculations = [
        vmc(
            walkers     =   1,
            warmupsteps =  20,
            blocks      = 200,
            steps       =  10,
            substeps    =   2,
            timestep    =  .4
            )
        ],
    dependencies = (conv,'other'),
    )

run_project()
