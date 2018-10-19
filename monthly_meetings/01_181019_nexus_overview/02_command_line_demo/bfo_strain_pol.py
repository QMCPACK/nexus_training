#! /usr/bin/env python

# imports from nexus and numpy
from numpy import arange,loadtxt
from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_vasp
from nexus import graph_sims

# nexus settings
settings(
    pseudo_dir    = './pseudopotentials',
    results       = '',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'oic5'
    )

# vasp job information
modules = '''
module ()
{
    eval `/opt/modules/3.1.6/bin/modulecmd bash $*`
}
module purge
module load composerxe/2011.10.319
module load mpi/openmpi-1.4.5-intel
'''
vasp = '/home/j1k/apps/vasp/vasp53/5.3.3/vasp.oic_opt'
vjob = job(nodes=2,hours=2,presub=modules,user_env=False,app=vasp)

# vasp inputs used by all runs
shared_inputs = obj(
    encut    = 550,
    ediff    = 1e-8,
    algo     = 'Normal',
    nelm     = 400,
    nwrite   = 2,
    prec     = 'Accurate',
    ismear   = -5,
    lwave    = False,
    lcharg   = False,
    lreal    = '.FALSE.',
    nbands   = 200,
    ispin    = 2,
    magmom   = 24*[0]+[-6,6,6,-6,6,-6,-6,6]+ 8*[0],
    ldau     = True,
    ldautype = 2,
    ldaul    = [-1, 2,-1],
    ldauu    = [ 0, 2, 0],
    ldauj    = [ 0, 0, 0],
    lmaxmix  = 4,
    lorbit   = 11,
    kpar     = 8,
    ncore    = 4,
    pseudos  = ['Bi.d.lda.POTCAR', # TITEL= PAW Bi_d 09Feb1998
                'Fe.pv.lda.POTCAR',# TITEL= PAW Fe_pv 03Mar1998
                'O.lda.POTCAR']    # TITEL= PAW O 31May2000
    )

# series of BFO cells to run
strains = arange(0.00,0.015,.005)  # strain %
shears  = arange(1.120,1.130,.005) # monoclinic shearing

# run geometry optimization, scf, and polarization 
# for all BFO cells 
for strain in strains:
    for shear in shears:
        path = 'st_{0:4.3f}_sh_{1:4.3f}/'.format(strain,shear)

        # strain and shearing cell parameters
        a = 7.80*(1.-strain)
        c = shear*a
        k = 0.0348994967*c
    
        # physical system for strained cell
        system = generate_physical_system(
            axes   = [[a,0,0],
                      [0,a,0],
                      [k,0,c]],
            elem   = 24*['O']+8*['Fe']+8*['Bi'],
            posu   = loadtxt('BFO_upos.dat'),
            units  = 'A',
            kgrid  = (2,2,2),
            kshift = (0,0,0)
            )

        # geometry relaxation run
        geo = generate_vasp(
            # standard inputs
            identifier   = 'geo',
            path         = path,
            job          = vjob,
            system       = system,
            # vasp inputs
            ibrion       = 2,
            isif         = 2,
            ediffg       = 1e-6,
            nsw          = 60,
            **shared_inputs
            )

        # scf run
        scf = generate_vasp(
            # standard inputs
            identifier   = 'scf',
            path         = path,
            job          = vjob,
            system       = system,
            dependencies = (geo,'other'),
            # vasp inputs
            **shared_inputs
            )

        # polarization run
        pol = generate_vasp(
            # standard inputs
            identifier   = 'pol',
            path         = path,
            job          = vjob,
            system       = system,
            dependencies = (scf,'other'),
            # vasp inputs
            lcalcpol     = True, 
            lepsilon     = True,
            **shared_inputs
            )
    #end for
#end for

run_project() 
