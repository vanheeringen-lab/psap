#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages
import versioneer

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("requirements.txt") as f:
    requirements = f.read().splitlines()


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
