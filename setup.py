#!/usr/bin/env python
from setuptools import find_packages, setup

setup_requires = [
    "cffi>=1.0.0",
]

install_requires = [
    "ase",
    "tqdm",
    "numpy",
    "h5py",
    "cffi",
    "ray",
]

setup(
    name="gmp-featurizer",
    version="0.1",
    description="Feature calculator for GMP features",
    author="Xiangyun Lei",
    author_email="ray.lei@tri.global",
    url="https://github.com/TRI-AMDD/GMP-featurizer",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    package_data={"": ["*.cpp", "*.h"]},
    python_requires=">=3.6, <4",
    setup_requires=setup_requires,
    install_requires=install_requires,
    cffi_modules=[
        "GMPFeaturizer/GMP/libgmp_builder.py:ffibuilder",
    ],
)
