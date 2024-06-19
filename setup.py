#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages
import codecs
import os.path

# Functions to pull the package version from init.py
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = ['pytest>=3', ]

setup(
    author="Tom Ellis",
    author_email='thomas.ellis@gmi.oeaw.ac.at',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Python tools for the epiclines ERC project",
    entry_points={
        'console_scripts': [
            'methlab=methlab.cli:main',
        ],
    },
    install_requires=['pandas', 'numpy','plotly','scikit-allel'],
    license="MIT license",
    #long_description=readme + '\n\n' + history,
    keywords='methlab',
    name='methlab',
    packages=find_packages(include=['methlab', 'methlab.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/ellisztamas/methlab',
    zip_safe=False,
    include_package_data=True,
    version=get_version("methlab/__init__.py"),
    package_data={
        "methlab": ["data/*.csv"]
    }
)
