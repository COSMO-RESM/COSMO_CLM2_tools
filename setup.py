import os
from setuptools import setup

def get_version():
    with open('COSMO_CLM2_tools/__init__.py') as f:
        for line in f:
            if line.startswith('__version__'):
                _, _, version = line.replace("'", '').split()
                break
    return version

setup(name='COSMO_CLM2_tools',
      version=get_version(),
      description="python based tools to set up a COSMO_CLM2 case",
      author="Matthieu Leclair",
      author_email="matthieu.leclair@env.ethz.ch",
      url="https://github.com/COSMO-RESM/COSMO-CLM2_tools",
      packages=['COSMO_CLM2_tools'],
      entry_points={'console_scripts': ['cc2_create_case = COSMO_CLM2_tools.create_case:create_case',
                                        'cc2_control_case = COSMO_CLM2_tools.control_case:control_case',
                                        'cc2_compile_clm = COSMO_CLM2_tools.compile_clm:compile_clm']},
      install_requires=['f90nml>=1.0.2']
)
