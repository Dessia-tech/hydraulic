
import os

# for some akward reason, put scripts before unittests tests  
scripts = ['2D/u.py',
           '3D/complex.py',
           '3D/hexagons.py',
           'thermal/complex.py',
           'thermal/unidimensional_medium.py'
           ]

for script_name in scripts:
    print('\n## Executing script {}'.format(script_name))
    exec(open(script_name).read())
