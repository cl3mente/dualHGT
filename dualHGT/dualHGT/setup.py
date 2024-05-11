from setuptools import setup, find_packages
setup(
name='dualHGT',
version='0.1.0',
author='Clemente Calabrese',
author_email='clemente.calabrese@gmail.com',
description='HGT inference based on third party software',
packages=find_packages(),
classifiers=[
'Programming Language :: Python :: 3',
'License :: OSI Approved :: MIT License',
'Operating System :: OS Independent',
],
python_requires='>=3.6',
)