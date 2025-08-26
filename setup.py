"""
    Bolt installation script
    Install bolt with `pip install -e .`, assuming you are in the top directory
"""

import subprocess
from setuptools import setup, find_packages

# TODO: try to make bolt install without editable

def build_lib():
    subprocess.check_call(["make"], cwd="bolt/")

build_lib()

setup(
    name='bolt',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[],
    entry_points={
        'console_scripts': [],
    },
    author='João Rebouças',
    author_email='joao.reboucas1@gmail.com',
    description='A library for linear cosmological perturbations',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/joaoreboucas1/bolt',
    classifiers=[
        'Development Status :: 1 - Planning',
        'Programming Language :: Python :: 3',
        'Programming Language :: C',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
    ],
    python_requires='>=3.6',
)