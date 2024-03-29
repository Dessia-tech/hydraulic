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

def version_from_git_describe(version):
    if version[0]=='v':
            version = version[1:]
            
    # PEP 440 compatibility
    number_commits_ahead = 0
    if '-' in version:
        version, number_commits_ahead, commit_hash = version.split('-')
        number_commits_ahead = int(number_commits_ahead)

    split_versions = version.split('.')
    if len(split_versions) == 4:
        suffix = split_versions[3]
        split_versions = split_versions[:3]
    else:
        suffix = None
    
    for suffix2 in ['a', 'b', 'rc']:
        if suffix2 in split_versions[-1]:
            if number_commits_ahead > 0:
                split_versions[-1] = str(split_versions[-1].split(suffix2)[0])
                split_versions[-1] = str(int(split_versions[-1])+1)
                future_version = '.'.join(split_versions)
                return '{}.dev{}'.format(future_version, number_commits_ahead)

            else:
                return '.'.join(split_versions)
    if number_commits_ahead > 0:
        split_versions[-1] = str(int(split_versions[-1])+1)
        split_versions = '.'.join(split_versions)
        return '{}.dev{}'.format(split_versions, number_commits_ahead)
    else:
        if suffix is not None:
            split_versions.append(suffix)

        return '.'.join(split_versions)
    
# Just testing if get_version works well
for v in ['v0.1.7.post2', 'v0.0.1-25-gaf0bf53', 'v0.0.1a2-25-gaf0bf53']:
    version_from_git_describe(v)
    
def get_version():
    # Return the version if it has been injected into the file by git-archive
    version = tag_re.search('$Format:%D$')
    if version:
        return version.group(1)

    d = dirname(__file__)
    
    if isdir(join(d, '.git')):
        cmd = 'git describe --tags'
        try:
            version = check_output(cmd.split()).decode().strip()[:]
            
        except CalledProcessError:
            raise RuntimeError('Unable to get version number from git tags')
        
        return version_from_git_describe(version)
    else:
        # Extract the version from the PKG-INFO file.
        with open(join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)
            
    # print('version', version)
    return version


setup(name='hydraulic',
      version = get_version(),
#      setup_requires=['setuptools_scm'],
      description=' A 1D hydraulic simulator in python. ',
      long_description=readme(),
      long_description_content_type='text/markdown',
      keywords='hydraulics, hydraulic',
      url='https://github.com/Dessia_tech/hydraulic',
      author='Dessia Technologies',
      author_email='root@dessia.tech',
      license='LGPL v3',
      packages=['hydraulic'],#,'volmdlr.primitives2D','volmdlr.primitives3D','volmdlr.geometry'],
      package_dir={},
      install_requires=['numpy', 'matplotlib', 'scipy', 'dessia_common>=0.10.0', 'volmdlr>0.8.0'],
      classifiers=['Topic :: Scientific/Engineering','Development Status :: 3 - Alpha'],
      )
