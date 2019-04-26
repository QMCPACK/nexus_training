#! /usr/bin/env python

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_quantum_package

# note: you must source the QP config file before running this script
#   source /your/path/to/quantum_package.rc

settings(
    results       = '',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws12',
    qprc          = \
'/home/j1k/apps/quantum_package/qp2-2.0.0-beta/quantum_package.rc',
    )

scf_job = job(cores=12,threads=12)

system = generate_physical_system(
    structure = 'H2O.xyz',
    )

scf = generate_quantum_package(
    identifier   = 'hf',        # log output goes to hf.out
    path         = 'h2o_ae_hf', # directory to run in
    job          = scf_job,
    system       = system,
    prefix       = 'h2o',       # create/use h2o.ezfio
    run_type     = 'scf',       # qprun scf h2o.ezfio
    ao_basis     = 'cc-pvtz',   # use cc-pvtz basis
    )

run_project()
