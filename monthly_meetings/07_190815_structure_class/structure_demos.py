

import os
import numpy as np
import matplotlib.pyplot as plt

from structure import generate_structure, read_structure
from testing import object_eq,value_eq

plt.rcParams.update({'legend.fontsize':14,'figure.facecolor':'white',
                     'figure.subplot.hspace':0.,'axes.labelsize':16,
                     'xtick.labelsize':14,'ytick.labelsize':14})


#=====================================================================#
#                       Structure generation                          #
#=====================================================================#

def demo_gen_h2o():
    """
    Create an H2O molecule.
    """

    h2o = generate_structure(
        elem  = ['O','H','H'], 
        pos   = [[0.000000, 0.000000, 0.000000],
                 [0.000000,-0.757160, 0.586260],
                 [0.000000, 0.757160, 0.586260]],
        units = 'A', # Angstroms
        )

    # print important internal attributes
    print('\nImportant internal data:')
    print('\nelem:     # Atoms in the cell')
    print(h2o.elem)
    print('\npos:      # Atomic positions')
    print(h2o.pos)
    print('\nunits:    # Angstrom (A) or Bohr (B)')
    print(h2o.units)
#end def demo_h2o_gen



def demo_gen_diamond_direct():
    """
    Create a conventional cell of diamond directly.
    """

    diamond = generate_structure(
        units = 'A',
        axes  = [[3.57, 0.00, 0.00],
                 [0.00, 3.57, 0.00],
                 [0.00, 0.00, 3.57]],
        elem  = 8*['C'],
        posu  = [[0.00, 0.00, 0.00],
                 [0.25, 0.25, 0.25],
                 [0.00, 0.50, 0.50],
                 [0.25, 0.75, 0.75],
                 [0.50, 0.00, 0.50],
                 [0.75, 0.25, 0.75],
                 [0.50, 0.50, 0.00],
                 [0.75, 0.75, 0.25]],
        )

    print('\nImportant internal data:')
    print('\nunits:    # Angstrom (A) or Bohr (B)')
    print(diamond.units)
    print('\nbconds:   # Boundary conditions in each direction')
    print(diamond.bconds)
    print('\naxes:     # Axes of the cell')
    print(diamond.axes)
    print('\ncenter:   # Center of the cell')
    print(diamond.center)
    print('\nkaxes:    # Reciprocal space axes of the cell')
    print(diamond.kaxes)
    print('\nelem:     # Atoms in the cell')
    print(diamond.elem)
    print('\npos:      # Atomic positions')
    print(diamond.pos)
    print('\nkpoints:  # List of k-points')
    print(diamond.kpoints)
    print('\nkweights: # Weight of each k-point')
    print(diamond.kweights)
    print('\nfrozen:   # Atoms not subject to relaxation or dynamics')
    print(diamond.frozen)
    print('\nfolded_structure: # Structure before tiling')
    print(diamond.folded_structure)

#end def demo_gen_diamond_direct



def demo_gen_diamond_lattice():
    """
    Create a conventional cell of diamond from lattice information.
    """

    diamond = generate_structure(
        lattice   = 'cubic',        # cubic tetragonal orthorhombic rhombohedral
                                    # hexagonal triclinic monoclinic
        cell      = 'conventional', # primitive, conventional
        centering = 'F',            # P A B C I R F
        constants = 3.57,           # a,b,c,alpha,beta,gamma
        units     = 'A',            # A or B
        atoms     = 'C',            # species in primitive cell
        basis     = [[0,0,0],       # basis vectors (optional)
                     [.25,.25,.25]],
        )

    # print structure data
    print(diamond)
#end def demo_gen_diamond_lattice



def demo_gen_graphene():
    """
    Create a graphene cell using stored information.
    """

    graphene = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        )

    # print structure data
    print(graphene)
#end def demo_gen_graphene




#=====================================================================#
#                   Writing structures to files                       #
#=====================================================================#

def demo_write():
    """
    Write conventional diamond cell to XYZ, XSF, and POSCAR formats.
    """

    d8 = generate_structure(
        structure = 'diamond',
        cell      = 'conv',
        )

    # Write an XYZ file
    d8.write('diamond8.xyz') # loses cell info, appropriate for molecules
    
    # Write an XSF file
    d8.write('diamond8.xsf')

    # Write a POSCAR file
    d8.write('diamond8.POSCAR')

    # Print the contents of each file
    print('\nXYZ file:')
    os.system('cat diamond8.xyz')

    print('\nXSF file:')
    os.system('cat diamond8.xsf')

    print('\nPOSCAR file:')
    os.system('cat diamond8.POSCAR')

#end def demo_write




#=====================================================================#
#                  Reading structures from files                      #
#=====================================================================#

def demo_read():
    """
    Read conventional diamond cell from XYZ, XSF, and POSCAR formats.
    """

    # Read an XYZ file
    d8_xyz = read_structure('diamond8.xyz') # no cell info
    
    # Read an XSF file
    d8_xsf = read_structure('diamond8.xsf')

    # Read a POSCAR file
    d8_poscar = read_structure('diamond8.POSCAR')

    # Compare with the reference structure
    d8 = generate_structure(
        structure = 'diamond',
        cell      = 'conv',
        )

    def struct_same(s1,s2,axes=True):
        same = True
        same &= value_eq(s1.elem,s2.elem)
        same &= value_eq(s1.pos,s2.pos)
        if axes:
            same &= value_eq(s1.axes,s2.axes)
        #end if
        return same
    #end def struct_same

    print('\nXYZ valid? {}'.format(struct_same(d8_xyz,d8,axes=False)))

    print('\nXSF valid? {}'.format(struct_same(d8_xsf,d8)))

    print('\nPOSCAR valid? {}'.format(struct_same(d8_poscar,d8)))

#end def demo_read



def demo_read_cif():
    """
    Read La2CuO4 structure from a CIF file.
    """

    # Note: this demo requires PyCifRW

    # Read from CIF file
    s = read_structure('La2CuO4_ICSD69312.cif')

    # Print out the structure
    print('La2CuO4 in POSCAR format:')
    print(s.write_poscar())

#end def demo_read_cif




#=====================================================================#
#                    Making molecular systems                         #
#=====================================================================#

def demo_water_box():
    """
    Make simulation cells for water molecule and dimer.
    """

    h2o = generate_structure(
        elem  = ['O','H','H'], 
        pos   = [[0.000000, 0.000000, 0.000000],
                 [0.000000,-0.757160, 0.586260],
                 [0.000000, 0.757160, 0.586260]],
        units = 'A', # Angstroms
        )

    # add a box by hand to the water molecule
    h2o_diy = h2o.copy()
    h2o_diy.set_axes(8.0*np.eye(3))

    # automatically add a bounding box to the water molecule
    h2o_auto = h2o.copy()
    h2o_auto.bounding_box(box='cubic',minsize=8.0)

    # make a water dimer configuration
    h2o1 = h2o.copy()
    h2o1.bounding_box(box='cubic',minsize=8.0)
    h2o2 = h2o1.copy()
    h2o1.translate((0,0,-1.1))
    h2o2.translate((0,0, 1.1))
    h2o_dimer = h2o1.copy()
    h2o_dimer.incorporate(h2o2)

    # plot the water molecule with DIY box
    plt.figure(tight_layout=True)
    h2o_diy.plot2d_ax(1,2,'k',lw=2)
    h2o_diy.plot2d_pos(1,2,'bo')
    plt.axis('equal')
    plt.xlabel('y (A)')
    plt.ylabel('z (A)')
    plt.title('demo water box: DIY box')

    # plot the water molecule with an automatic box
    plt.figure(tight_layout=True)
    h2o_auto.plot2d_ax(1,2,'k',lw=2)
    h2o_auto.plot2d_pos(1,2,'bo')
    plt.axis('equal')
    plt.xlabel('y (A)')
    plt.ylabel('z (A)')
    plt.title('demo water box: automatic box')

    # plot the water dimer
    plt.figure(tight_layout=True)
    h2o_dimer.plot2d_ax(1,2,'k',lw=2)
    h2o_dimer.plot2d_pos(1,2,'bo')
    plt.axis('equal')
    plt.xlabel('y (A)')
    plt.ylabel('z (A)')
    plt.title('demo water box: water dimer')

    plt.show()
#end def demo_water_box



def demo_coronene_box():
    """
    Place coronene in a well apportioned box.
    """

    s = read_structure('./coronene.xyz')

    # make a bounding box that is at least 5 A from the nearest atom
    s.bounding_box(mindist=5.0)

    # print the resulting box
    print(s.axes)

    # plot the coronene molecule with an automatic box
    plt.figure(tight_layout=True)
    s.plot2d_ax(0,1,'k',lw=2)
    s.plot2d_pos(0,1,'bo')
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo coronene box: box with minimum distance')
    plt.show()
#end def demo_coronene_box




#=====================================================================#
#                          Tiling structures                          #
#=====================================================================#

def demo_tiling():
    """
    Make a 4x4 tiled cell of graphene and plot it.
    """

    g11 = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        )

    # Get the tiled supercell
    print('Tiling the primitive cell')
    g44 = g11.tile(4,4,1)

    # Tiled cells contain the originating ("folded") ones
    g11p = g44.folded_structure

    # compare positions of the 1x1 cells
    folded_same = g11.distances(g11p).max()<1e-6
    print('Folded structure agrees with original? {}'.format(folded_same))
    
    # By default, folded structures are used in DFT and tiled ones in QMC.
    # You can discard the folded structure this way:
    print('\nRemoving folded structure')
    g44.remove_folded()
    print('Folded structure is gone? {}'.format(not g44.has_folded()))

    # Plot the tiled structure
    plt.figure(tight_layout=True)
    g44.plot2d_ax(0,1,'k',lw=2)
    g11.plot2d_ax(0,1,'r',lw=2)
    g44.plot2d_pos(0,1,'bo')
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo tiling: 4x4 tiled graphene')
    plt.show()
#end def demo_tiling



def demo_matrix_tiling():
    """
    Make a matrix tiled cell of graphene and plot it.
    """

    g11 = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        )

    tmat = [[3, 0, 0],
            [3, 6, 0],
            [0, 0, 1]]

    gmat = g11.tile(tmat)

    # For any two structures related by tiling, tiling matrix can be found:
    print('\nTiling matrix between primitive and supercell:')
    print(gmat.tilematrix(g11))

    # Plot the tiled structure
    plt.figure(tight_layout=True)
    gmat.plot2d_ax(0,1,'k',lw=2)
    g11.plot2d_ax(0,1,'r',lw=2)
    gmat.plot2d_pos(0,1,'bo')
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo tiling: matrix tiled graphene')
    plt.show()
#end def demo_matrix_tiling



def demo_opt_tiling():
    """
    Search for the "best" tiling matrix for a graphene supercell.
    """

    g11 = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        )

    # Look for an optimal cell 20x larger than the primitive cell
    gopt = g11.tile_opt(20)

    print('\nOptimal(?) tiling matrix between primitive and supercell:')
    print(gopt.tilematrix(g11))

    # Plot the tiled structure
    plt.figure(tight_layout=True)
    gopt.plot2d_ax(0,1,'k',lw=2)
    g11.plot2d_ax(0,1,'r',lw=2)
    gopt.plot2d_pos(0,1,'bo')
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo tiling: optimally tiled graphene')
    plt.show()
#end def demo_opt_tiling




#=====================================================================#
#                     Finding primitive cells                         #
#=====================================================================#

def demo_primitive_search():
    """
    Find the primitive cell given a supercell.
    """

    # Note: this demo requires SeeKpath

    d2 = generate_structure(
        structure = 'diamond',
        cell      = 'prim',
        )

    # Construct the cubic 64 atom supercell
    d64 = d2.tile_opt(32)

    # Remove all traces of the 2 atom cell, supercell is all that remains
    d64.remove_folded()
    del d2

    # Find the primitive cell from the supercell
    dprim = d64.primitive()

    # Print the supercell axes
    print('\nSupercell axes:\n{}'.format(d64.axes))

    # Print the primitive cell axes
    print('\nPrimitive cell axes:\n{}'.format(dprim.axes))

    # Print the tiling matrix
    print('\nTiling matrix:\n{}'.format(d64.tilematrix(dprim)))

#end def demo_primitive_search




#=====================================================================#
#                          Adding k-points                            #
#=====================================================================#

def demo_mp_kpoints():
    """
    Add k-points according to Monkhorst-Pack grid.
    """

    g11 = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        )

    # Perform a 4x4 tiling of the primitive cell
    g44 = g11.tile(4,4,1)

    # Add a Gamma-centered 2x2 Monkhorst-Pack grid
    g44g = g44.copy()
    g11g = g44g.folded_structure

    g44g.add_kmesh(kgrid=(2,2,1),kshift=(0,0,0))

    # Add a shifted 2x2 Monkhorst-Pack grid
    g44s = g44.copy()
    g11s = g44s.folded_structure

    g44s.add_kmesh(kgrid=(2,2,1),kshift=(0.5,0.5,0))

    # Get the mapping between supercell and primitive cell k-points
    print('Monk-horst pack grid')
    kmap = g44s.kmap()
    print('# of k-points in supercell: {}'.format(len(g44s.kpoints)))
    print('# of k-points in primcell : {}'.format(len(g11s.kpoints)))
    print('index mapping from super to prim:\n{}'.format(str(kmap).replace('  ','')))

    # Plot kpoints for the 2x2 unshifted MP mesh
    plt.figure(tight_layout=True)
    g11g.plot2d_kax(0,1,'r',lw=2,label='primcell')
    g44g.plot2d_kax(0,1,'k',lw=2,label='supercell')
    g11g.plot2d_kp(0,1,'g.')
    g44g.plot2d_kp(0,1,'bo')
    plt.legend()
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo kpoints: unshifted 2x2 MP mesh in 4x4 graphene')

    # Plot kpoints for the 2x2 shifted MP mesh
    plt.figure(tight_layout=True)
    g11s.plot2d_kax(0,1,'r',lw=2,label='primcell')
    g44s.plot2d_kax(0,1,'k',lw=2,label='supercell')
    g11s.plot2d_kp(0,1,'g.')
    g44s.plot2d_kp(0,1,'bo')
    plt.legend()
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo kpoints: shifted 2x2 MP mesh in 4x4 graphene')

    plt.show()
#end def demo_mp_kpoints



def demo_random_kpoints():
    """
    Add random k-points.
    """

    g11 = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        )

    # Perform a 4x4 tiling of the primitive cell
    g44 = g11.tile(4,4,1)

    # Add random kpoints instead of an MP mesh
    rand = np.random.uniform(size=(20,2))
    kpoints = np.dot(rand,g44.kaxes[0:2])
    g44.add_kpoints(kpoints)
    g11 = g44.folded_structure

    # Plot the random kpoints
    plt.figure(tight_layout=True)
    g11.plot2d_kax(0,1,'r',lw=2,label='primcell')
    g44.plot2d_kax(0,1,'k',lw=2,label='supercell')
    g11.plot2d_kp(0,1,'g.')
    g44.plot2d_kp(0,1,'bo')
    plt.legend()
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo kpoints: random kpoints in 4x4 graphene')

    plt.show()
#end def demo_random_kpoints



def demo_symm_kpoints():
    """
    Add symmetrized Monkhorst-Pack kpoints.
    """

    # Note: this demo requires spglib

    g44 = generate_structure(
        structure  = 'graphene',
        cell       = 'prim',
        tiling     = (4,4,1),
        kgrid      = (4,4,1),
        kshift     = (0,0,0),
        symm_kgrid = True,
        )

    g11 = g44.folded_structure

    # Print the symmtrized kpoints
    print('Symmetrized k-points:\n{}'.format(g44.kpoints_unit()))
    print('Symmetrized weights:\n{}'.format(g44.kweights))

    # Plot the symmetrized kpoints
    plt.figure(tight_layout=True)
    g11.plot2d_kax(0,1,'r',lw=2,label='primcell')
    g44.plot2d_kax(0,1,'k',lw=2,label='supercell')
    g11.plot2d_kp(0,1,'g.')
    g44.plot2d_kp(0,1,'bo')
    plt.legend()
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo kpoints: symmetrized kpoints in 4x4 graphene')

    plt.show()

#end def demo_symm_kpoints




#=====================================================================#
#                      Getting cell information                       #
#=====================================================================#

def demo_rwigner():
    """
    Show Wigner and inscribing radii graphically.
    """

    s = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        )

    # get the tiled supercell
    st = s.tile(4,4,1)

    # get the wigner radius
    rw = st.rwigner()

    # get the incribing radius
    ri = st.rinscribe()

    # Print the radii
    print('\nWigner radius: {}'.format(rw))
    print('\nInscribing radius: {}'.format(ri))
    
    # plot the tiled structure
    plt.figure(tight_layout=True)
    st.plot2d_ax(0,1,'k',lw=2)
    st.plot2d_pos(0,1,'bo')
    t = np.linspace(0,2*np.pi,400,endpoint=True)
    sin = np.sin(t)
    cos = np.cos(t)
    xc,yc,zc = st.center
    plt.plot(rw*cos+xc,rw*sin+yc,'r-',lw=2,label='r wigner')
    plt.plot(ri*cos+xc,ri*sin+yc,'g-',lw=2,label='r inscribe')
    plt.axis('equal')
    plt.legend(loc='upper right')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo rwigner and rinscribe')
    plt.show()
#end def demo_rwigner



def demo_cell_constants():
    """
    Compute the cell volume, Madelung constant, and Makov Payne correction
    """

    d2 = generate_structure(
        structure = 'diamond',
        cell      = 'prim',
        )

    # Construct the cubic 64 atom supercell
    d64 = d2.tile_opt(32)

    # Print the volume
    print('\nVolume of 64 atom cell (A): {}'.format(d64.volume()))
    
    # Print the Madelung constant
    #   ( see equation 7 in PRB 78 125106 (2008) )
    print('\nMadelung constant (Ha): {}'.format(d64.madelung()))

    # Print Makov-Payne corrections
    mp_q1 = d64.makov_payne(q=1,eps=5.68)
    mp_q2 = d64.makov_payne(q=2,eps=5.68)
    print('\nMakov-Payne corr. (Ha) for charged defect (q=1): {}'.format(mp_q1))
    print('\nMakov-Payne corr. (Ha) for charged defect (q=2): {}'.format(mp_q2))

#end def demo_cell_constants



def demo_unit_coords():
    """
    Get atomic positions and k-points in unit coordinates.
    """

    s = generate_structure(
        structure = 'diamond',
        cell      = 'prim',
        tiling    = (2,2,2),
        kgrid     = (2,2,2),
        kshift    = (0.5,0.5,0.5),
        )

    print('\nAtomic positions in unit coordinates:\n{}'.format(s.pos_unit()))

    print('\nk-points in unit coordinates:\n{}'.format(s.kpoints_unit()))

#end def demo_unit_coords




#=====================================================================#
#      Calculating minimum image distances and nearest neighbors      #
#=====================================================================#

def demo_min_image_distances():
    """
    Compute minimum image distances between nearest neighbors.
    """

    g = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (4,4,1),
        )

    # Get the neighbor (index) table, along with sorted distance 
    # and displacement tables.
    nt,dt,vt = g.neighbor_table(distances=True,vectors=True)

    # Restrict to 3 nearest neighbors (not including self at index 0)
    nt = nt[:,1:4]
    dt = dt[:,1:4]
    vt = vt[:,1:4]
    
    # Print the nearest neighbor distance table
    print('\nDistances to 3 nearest neighbors for each atom:\n{}'.format(dt))
    
    # Print the nearest neighbor indices
    print('\nIndices of 3 nearest neighbors for each atom:\n{}'.format(nt))

    # Plot the displacement vectors
    plt.figure(tight_layout=True)
    g.plot2d_ax(0,1,'k',lw=2)
    for r,v in zip(g.pos,vt):
        for d in v:
            plt.plot([r[0],r[0]+d[0]],[r[1],r[1]+d[1]],'g-')
        #end for
    #end for
    g.plot2d_pos(0,1,'bo')
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo min image dist: distances between neighboring atoms')
    plt.show()
#end def demo_min_image_distances



def demo_compare_structures():
    """
    Compute distances and displacements between two structures.
    """

    g = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (4,4,1),
        )

    # make random displacements to represent a second stucture
    gr = g.copy()
    gr.pos += 1.5*(np.random.uniform(size=gr.pos.shape)-0.5)
    # put the atoms back in the box
    gr.recenter()

    # get distance information, but not over all pairs
    dist,disp = g.min_image_distances(gr,vectors=True,pairs=False)

    # print out the distance results
    print('\nPairwise distances between the structures:\n{}'.format(dist))

    # Plot the displacement vectors
    plt.figure(tight_layout=True)
    g.plot2d_ax(0,1,'k',lw=2)
    for r,d in zip(gr.pos,disp):
        plt.plot([r[0],r[0]+d[0]],[r[1],r[1]+d[1]],'g-')
    #end for
    g.plot2d_pos(0,1,'bo')
    gr.plot2d_pos(0,1,'r.')
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo compare structures: distances between neighboring atoms')
    plt.show()

#end def demo_compare_structures





#=====================================================================#
#            Freezing sets of atoms to prevent relaxation             #
#=====================================================================#

def demo_freeze():
    """
    Freeze sets of atoms to prevent relaxation.
    """

    g = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (6,6,1),
        )
    g.recenter((0,0,0))

    # Select atoms to be frozen
    outer_atoms = np.sqrt((g.pos**2).sum(1)) > 3.2

    # Freeze the atoms
    g.freeze(outer_atoms)

    # Get a mask array identifying the frozen atoms
    frozen = g.is_frozen()

    # Plot the frozen and movable atoms
    plt.figure(tight_layout=True)
    g.plot2d_ax(0,1,'k',lw=2)
    t = 2*np.pi*np.linspace(0,1,400)
    plt.plot(3.2*np.cos(t),3.2*np.sin(t),'g--')
    plt.plot(g.pos[frozen,0],g.pos[frozen,1],'bo',label='frozen')
    plt.plot(g.pos[~frozen,0],g.pos[~frozen,1],'ro',label='movable')
    plt.axis('equal')
    plt.legend()
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo freeze: prevent atoms from relaxing')
    plt.show()
#end def demo_freeze




#=====================================================================#
#       Embedding relaxed structures in larger pristine cells         #
#=====================================================================#

def demo_embed():
    """
    Embed a "relaxed" structure in a larger pristine cell.
    """

    center = (0,0,0)

    g = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (4,4,1),
        )
    g.recenter(center)

    # Represent the "relaxed" cell
    gr = g.copy()
    npos = len(gr.pos)
    dr = gr.min_image_vectors(center)
    dr.shape = npos,3
    r = np.linalg.norm(dr,axis=1)
    dilation = 2*r*np.exp(-r)
    for i in range(npos):
        gr.pos[i] += dilation[i]/r[i]*dr[i]
    #end for

    # Represent the unrelaxed large cell
    gl = generate_structure(
        structure = 'graphene',
        cell      = 'rect',
        tiling    = (8,4,1),
        )
    gl.recenter(center)

    # Embed the relaxed cell in the large unrelaxed cell
    ge = gl.copy()
    ge.embed(gr)

    # Plot the relaxed structure
    plt.figure(tight_layout=True)
    gr.plot2d_ax(0,1,'k',lw=2)
    g.plot2d_pos(0,1,'c.')
    gr.plot2d_pos(0,1,'bo')
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo embed: small relaxed in large pristine')

    # Plot the large unrelaxed cell
    plt.figure(tight_layout=True)
    gl.plot2d_ax(0,1,'k',lw=2)
    gl.plot2d_pos(0,1,'bo')
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo embed: small relaxed in large pristine')

    # Plot the large unrelaxed cell
    plt.figure(tight_layout=True)
    gl.plot2d_ax(0,1,'k',lw=2)
    gr.plot2d_ax(0,1,'g',lw=2)
    gr.plot2d_pos(0,1,'ro')
    gl.plot2d_pos(0,1,'b.')
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo embed: small relaxed in large pristine')

    # Plot the large unrelaxed cell
    plt.figure(tight_layout=True)
    ge.plot2d_ax(0,1,'k',lw=2)
    ge.plot2d_pos(0,1,'bo')
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo embed: small relaxed in large pristine')

    plt.show()
#end def demo_embed




#=====================================================================#
#               Interpolation for nudged elastic band                 #
#=====================================================================#

def demo_interpolate():
    """
    Interpolate between two "relaxed" structures for NEB initialization.
    """
    
    g = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (4,4,1),
        )
    g11 = g.folded_structure

    # Make chromium atom positions (pretend like they are above the surface)
    npos = np.dot((1./3,2./3,0),g11.axes)
    npos1 = npos+3*g11.axes[0]
    npos2 = npos1+g11.axes[0]+g11.axes[1]

    # "Relaxed" structure with additional atom on one ring
    gr1 = g.copy()
    gr1.add_atoms(['Cr'],[npos1])

    # "Relaxed" structure with additional atom on neighboring ring
    gr2 = g.copy()
    gr2.add_atoms(['Cr'],[npos2])
    gr2.recenter()

    # Interpolate between the two structures within min. image convention
    spath = gr1.interpolate(gr2,10)
    
    # Plot the (linear) interpolation path
    plt.figure(tight_layout=True)
    for gr in spath:
        gr.plot2d_ax(0,1,'k',lw=2)
        gr.plot2d_pos(0,1,'bo')
        npos = gr.pos[-1]
        plt.plot([npos[0]],[npos[1]],'go')
    #end for
    plt.axis('equal')
    plt.xlabel('x (A)')
    plt.ylabel('y (A)')
    plt.title('demo interpolate: paths for NEB')

    plt.show()
#end def demo_interpolate






from developer import error
import inspect
def demo(name):
    """
    Function used to run the demos.
    """
    name = 'demo_'+name
    g = globals()
    if name not in g:
        error('requested "{}" does not exist'.format(name))
    #end if
    demo_function = g[name]
    def write_header(s):
        print('')
        print(s)
        print(len(s)*'-')
    #end def write_header
    print(80*'=')
    write_header('demo code:')
    print(inspect.getsource(demo_function))
    write_header('demo execution:')
    demo_function()
#end def demo
