try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='COSMO_CLM2_tools',
      version='0.1',
      description="python based tools to set up a cOSMO_CLM2 case",
      author="Matthieu Leclair",
      author_email="matthieu.leclair@env.ethz.ch",
      url="https://github.com/COSMO-RESM/COSMO-CLM2_tools",
      packages=['COSMO_CLM2_tools'],
      entry_points={'console_scripts': ['COSMO_CLM2_tools.cc2_create_case = cosmo_clm2:create_new_case']},
)
