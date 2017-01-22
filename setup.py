from setuptools import setup, find_packages

setup(name='tanner',
      version='0.1',
      description='Analysis code related to the tanner project',
      url='http://github.com/Schork-Lab/tanner',
      author='Kunal Bhutani',
      author_email='kbhutani@ucsd.edu',
      license='MIT',
      packages=find_packages(),
      zip_safe=False,
      entry_points={'console_scripts': ['metabolome-fit=tanner.analysis.metabolome:bayesian_fit'],})
