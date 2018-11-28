#! /usr/bin/env python

from nexus import settings,job,run_project
from nexus import generate_structure
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack,vmc
from nexus import bundle

settings(
    pseudo_dir    = './pseudopotentials',
    results       = '',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws16',
    #machine       = 'eos',
    #account       = 'mat151',
    )

on_desktop = settings.machine=='ws16'
on_cluster = settings.machine=='eos'

if on_desktop:
    scf_job  = job(cores=16,app='pw.x')
    conv_job = job(cores=1,app='pw2qmcpack.x')
    qmc_job  = job(cores=1) # fake/placeholder job
elif on_cluster:
    scf_job  = job(nodes=1) # fake/placeholder job
    conv_job = job(nodes=1) # fake/placeholder job
    modules = '''
source $MODULESHOME/init/bash
if (echo $LOADEDMODULES | grep -q pgi)
then
module unload PrgEnv-pgi
fi
if (echo $LOADEDMODULES | grep -q gnu)
then
module unload PrgEnv-gnu
fi
if (echo $LOADEDMODULES | grep -q hdf5)
then
module unload cray-hdf5
fi
if (echo $LOADEDMODULES | grep -q libsci)
then
module unload cray-libsci
fi
module load PrgEnv-intel
module unload gcc
module load gcc
module load cray-hdf5-parallel
module load fftw
module load boost
module load subversion
module load cmake
'''
    qmcpack = '/ccs/home/jtkrogel/apps/eos/qmcpack/qmcpack-3.5.0/qmcpack_soa_real'
    qmc_job  = job(nodes   = 1,
                   threads = 4,
                   minutes = 30,
                   queue   = 'debug',
                   presub  = modules,
                   app     = qmcpack)
else:
    print 'not on desktop or cluster!'
    exit()
#end if


scalings = [0.95,1.00,1.05,1.10]

qmc_runs = []

for scaling in scalings:

    dia16_struct = generate_structure(
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
        )

    dia16_struct.rescale(scaling)

    dia16 = generate_physical_system(
        structure = dia16_struct,
        C         = 4
        )

    basepath = 'diamond/scale_{0:3.2f}'.format(scaling)

    scf = generate_pwscf(
        identifier   = 'scf',
        path         = basepath+'/scf',
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
        path         = basepath+'/scf',
        job          = conv_job,
        write_psir   = False,
        dependencies = (scf,'orbitals'),
        )

    qmc = generate_qmcpack(
        block        = not on_cluster,
        identifier   = 'vmc',
        path         = basepath+'/vmc',
        job          = qmc_job,
        input_type   = 'basic',
        system       = dia16,
        pseudos      = ['C.BFD.xml'],
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
        dependencies = (conv,'orbitals'),
        )

    qmc_runs.append(qmc)

#end for

if on_cluster:
    bqmc = bundle(qmc_runs)
#end if

run_project()
