import os
from setuptools import setup
# try:
#     from setuptools import setup
# except ImportError:
#     from distutils.core import setup
    
setup(name='COSMO_CLM2_tools',
      version='0.1',
      description="python based tools to set up a cOSMO_CLM2 case",
      author="Matthieu Leclair",
      author_email="matthieu.leclair@env.ethz.ch",
      url="https://github.com/COSMO-RESM/COSMO-CLM2_tools",
      packages=['COSMO_CLM2_tools'],
      entry_points={'console_scripts': ['cc2_create_case = COSMO_CLM2_tools.cosmo_clm2:create_new_case',
                                        'cc2_control_case = COSMO_CLM2_tools.cosmo_clm2:control_case'],},
      install_requires=['f90nml>=1.0.2']
)
