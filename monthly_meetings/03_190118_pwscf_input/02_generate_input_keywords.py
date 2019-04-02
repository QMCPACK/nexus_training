#! /usr/bin/env python

from nexus import generate_pwscf_input


print '''

pw = generate_pwscf_input(
    selector = 'generic',
    # what can go here?
    )
'''


print  '\n#most of the answer:'

from pwscf_input import PwscfInputBase
all_variables = PwscfInputBase.all_variables

def render(vlist):
    s = ''
    for v in sorted(vlist):
        if len(s+v)<80:
            s += '  '+v
        else:
            print s
            s = ''
        #end if
    #end for
    if len(s)>0:
        print s
    #end if
#end def render

print
print
print '==================================='
print 'PWSCF namelist variables supported:'
print '==================================='
render(all_variables)

print
print
print '=============================='
print 'Namelist variables by namelist'
print '(compare to INPUT_PW.html)'
print '=============================='
from pwscf_input import section_classes # namelist classes
for namelist in section_classes:
    print
    print namelist.__name__
    print '-------------------'
    render(namelist.variables)
#end for


print
print
print
print '================='
print 'Breakdown by type'
print '================='
print
print 'integers:'
print '---------'
render(PwscfInputBase.ints)
print
print 'floats:'
print '-------'
render(PwscfInputBase.floats)
print
print 'strings:'
print '--------'
render(PwscfInputBase.strs)
print
print 'booleans:'
print '---------'
render(PwscfInputBase.bools)
print
print 'arrays:'
print '------'
render(PwscfInputBase.real_arrays)


print
print
print
print '==============' 
print 'Valid examples'
print '==============' 
print '''
integers:
---------
generate_pwscf_input(
    nspin = 2,
    nbnd  = 63,
    )
'''
print '''
floats:
---------
generate_pwscf_input(
    ecutwfc     = 350.0,
    mixing_beta = 0.5,
    )
'''
print '''
strings:
--------
generate_pwscf_input(
    calculation = 'scf',
    input_dft   = 'lda',
    )
'''
print '''
booleans:
---------
generate_pwscf_input(
    lda_plus_u = True,
    wf_collect = False,
    )
'''
print '''
arrays:
-------
generate_pwscf_input(
    hubbard_u              = {1 : 3.1},
    starting_magnetization = {1 : 0.9},
    starting_ns_eigenvalue = {(1,2,1) : 0.0,
                              (2,2,1) : 0.0476060,
                              (3,2,1) : 0.0476060,
                              (4,2,1) : 0.9654373,
                              (5,2,1) : 0.9954307},
    )

generate_pwscf_input(
    hubbard_u              = obj(V1=3.5,V2=3.5),
    starting_magnetization = obj(V1=1.0,V2=-1.0),
    )
'''


print
print
print 
print 'This covers namelist input pretty well,'
print 'but what about "card" input?'
print
print 'More is needed, either from Nexus-specific'
print 'keyword inputs or data rich objects such'
print 'as PhysicalSystem.'
print
print 'Alternately one can read/compose input directly'
print '(bypass generate_pwscf_input).'

