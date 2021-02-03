#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()


setup(
    author="Tilman Schaefers",
    author_email='tilman.schaefers@ru.nl',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="CLI interface for the PSAP classifier, Mierlo, G. van. Predicting protein condensate formation using machine learning (Manuscript in Preparation).",
    entry_points={
        'console_scripts': [
            'psap=psap.cli:main',
        ],
    },
    license="MIT license",
    long_description=readme,
    keywords='psap',
    name='psap',
    packages=find_packages(include=['psap', 'psap.*']),
    install_requires=requirements,
    url='https://github.com/tilschaef/psap',
    version='0.1.3'
)
