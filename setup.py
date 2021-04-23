#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages
import versioneer

with open("README.rst") as readme_file:
    readme = readme_file.read()

requirements = [
    "pip==19.2.3",
    "wheel>=0.36.2",
    "scikit-learn>=-0.21.3",
    "twine>=1.14.0",
    "pandas>=1.0.1",
    "biopython>=1.73",
    "scipy>=1.2.0",
    "tqdm>=4.38.0",
    "seaborn>=0.11.1",
    "matplotlib>=3.3.4",
    "sklearn-json>=0.1.0",
    "versioneer>=0.19",
    "loguru>=0.5.3",
    "black>=20.8b1",
]



setup(
    author=["Juriaan Jansen", "Tilman Schaefers <tilman.schaefers@ru.nl>"],
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="CLI interface for the PSAP classifier. PSAP implements a RandomForest approach to predict the probability of proteins to mediate protein phase separation (PPS).",
    entry_points={
        "console_scripts": [
            "psap=psap.cli:main",
        ],
    },
    license="MIT license",
    long_description_content_type="text/x-rst",
    long_description=readme,
    keywords="psap",
    name="psap",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    install_requires=requirements,
    packages=find_packages(include=["psap", "psap.*"]),
    include_package_data=True,
    url="https://github.com/vanheeringen-lab/psap",
)
