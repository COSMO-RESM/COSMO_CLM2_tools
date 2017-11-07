from setuptools import setup

setup(name="COSMO_CLM2_tools",
      version="0.1",
      description="python based tools to set up a cOSMO_CLM2 case",
      author="Matthieu Leclair",
      author_email="matthieu.leclair@env.ethz.ch",
      url="https://git.iac.ethz.ch/leclairm/pv_cam_fv-dycore",
      entry_points={'console_scripts': ['cc2_create_case = cosmo_clm2:create_new_case']},
)
