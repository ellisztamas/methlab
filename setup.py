#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

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
            'align_filenames=epiclines_tools.cli:align',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    keywords='epiclines_tools',
    name='epiclines_tools',
    packages=find_packages(include=['epiclines_tools', 'epiclines_tools.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/ellisztamas/epiclines_tools',
    version='0.1.2',
    zip_safe=False,
    include_package_data=True,
    package_data={
        "epiclines_tools": ["data/*.csv"]
    }
)
