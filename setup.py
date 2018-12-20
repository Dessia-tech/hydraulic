# -*- coding: utf-8 -*-
"""
Setup install script for volmdlr

"""

from setuptools import setup
#from distutils.core import setup

from os.path import dirname, isdir, join
import re
from subprocess import CalledProcessError, check_output

tag_re = re.compile(r'\btag: %s([0-9][^,]*)\b')
version_re = re.compile('^Version: (.+)$', re.M)

def readme():
    with open('README.md') as f:
        return f.read()
    
def get_version():
    # Return the version if it has been injected into the file by git-archive
    version = tag_re.search('$Format:%D$')
    if version:
        return version.group(1)

    d = dirname(__file__)
    
    if isdir(join(d, '.git')):
        cmd = 'git describe --tags  --dirty'
        try:
            version = check_output(cmd.split()).decode().strip()[:]
        except CalledProcessError:
            version = 'v0.0.0'
#            raise RuntimeError('Unable to get version number from git tags')
        if version[0]=='v':
            version = version[1:]
#        print(version)
        # PEP 440 compatibility
        if '-' in version:
            if version.endswith('-dirty'):
                version = '.dev'.join(version.split('-')[:-1][:2])+'-dirty'
        else:
            version = '.dev'.join(version.split('-')[:2])

    else:
        # Extract the version from the PKG-INFO file.
        with open(join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)
            
#    # Writing to file
#    with open('powertransmission/version.py', 'w+') as vf:
#        vf.write("# -*- coding: utf-8 -*-\nversion = '{}'".format(version))
                 
    return version


setup(name='hydraulic',
      version = get_version(),
#      setup_requires=['setuptools_scm'],
      description=' A 1D hydraulic simulator in python. ',
      long_description=readme(),
      keywords='hydraulics, hydraulic',
      url='https://github.com/Dessia_tech/hydraulic',
      author='Dessia Technologies',
      author_email='root@dessia.tech',
      license='LGPL v3',
      packages=['hydraulic'],#,'volmdlr.primitives2D','volmdlr.primitives3D','volmdlr.geometry'],
      package_dir={},
      install_requires=['numpy', 'matplotlib', 'scipy', 'volmdlr'],
      classifiers=['Topic :: Scientific/Engineering','Development Status :: 3 - Alpha'],
      )
