from os import path
from setuptools import setup

from orthosnap.version import __version__

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

CLASSIFIERS = [
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Topic :: Scientific/Engineering",
]

REQUIRES = ["biopython>=1.81", "numpy>=1.24.0", "tqdm>=4.66.1"]

setup(
    name="orthosnap",
    description="orthosnap, identify orthologous subgroups of genes in large orthologous groups of genes.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Jacob L. Steenwyk",
    author_email="jlsteenwyk@gmail.com",
    url="https://github.com/jlsteenwyk/orthosnap",
    packages=["orthosnap"],
    classifiers=CLASSIFIERS,
    entry_points={"console_scripts": ["orthosnap = orthosnap.orthosnap:main"]},
    version=__version__,
    include_package_data=True,
    install_requires=REQUIRES,
)

# push new version to pypi
# rm -rf dist
# python3 setup.py sdist bdist_wheel --universal
# twine upload dist/* -r pypi
