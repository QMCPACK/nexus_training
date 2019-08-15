#!/usr/bin/env python

from pyscf.pbc import df, scf


### generated system text ###
from numpy import array
from pyscf.pbc import gto as gto_loc
cell = gto_loc.Cell()
cell.a             = '''
                     1.78500000   1.78500000   0.00000000
                     0.00000000   1.78500000   1.78500000
                     1.78500000   0.00000000   1.78500000
                     '''
cell.basis         = 'bfd-vdz'
cell.dimension     = 3
cell.ecp           = 'bfd'
cell.unit          = 'A'
cell.atom          = '''
                     C    0.00000000   0.00000000   0.00000000
                     C    0.89250000   0.89250000   0.89250000
                     '''
cell.drop_exponent = 0.1
cell.verbose       = 5
cell.charge        = 0
cell.spin          = 0
cell.build()
kpts = array([
    [0.0, 0.0, 0.0],
    [0.8799979421820149, 0.8799979421820149, -0.8799979421820149]])
### end generated system text ###



gdf = df.GDF(cell,kpts)
gdf.auxbasis = 'weigend'
gdf.build()

mf = scf.KRKS(cell,kpts).density_fit()
mf.xc      ='b3lyp'
mf.tol     = 1e-10
mf.exxdiv  = 'ewald'
mf.with_df = gdf
mf.kernel()

### generated conversion text ###
from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(cell,mf,'scf',kpts)
### end generated conversion text ###

