---
title: 'GMP-Featurizer: A parallelized Python package for efficiently computing the Gaussian Multipole features of atomic systems'
tags:
  - Python
  - C++
  - Parallelization
  - Machine Learning
  - Chemistry
  - Molecular Dynamics
authors:
  - name: Xiangyun Lei
    affiliation: 1
  - name: Joseph Montoya
    affiliation: 1
affiliations:
 - name: Toyota Research Institute, Los Altos, CA, USA
   index: 1
date: 4 April 2023
bibliography: paper.bib
---

# Summary

GMP-Featurizer is a lightweight, accurate, efficient, and scalable software package for calculating the Gaussian Multipole (GMP) features [@GMP] for a variety of atomic systems with elements across the periodic table. Starting from the GMP feature computation module from AmpTorch [@amptorch], the capability of GMP-Featurizer has since been greatly improved, including its accuracy and efficiency (please refer to the Overview section for details), as well as the ability to parallelize on different cores, even machines. Moreover, this Python package only has very few dependencies that are all standard python libraries, plus CFFI for C++ code interfacing and Ray [@Ray] for parallelization,  making it lightweight and robust. A set of unit tests are designed to ensure the reliability of its outputs. A set of extensive examples and tutorials, as well as two sets of pseudopotential files (needed for specifying the GMP feature set), are also included in this package for its users. Overall, this package is designed to serve as a standard implementation for chemical and material scientists who are interested in developing models based on GMP features. The source code for this package is freely available to the public under the Apache 2.0 license.

# Statement of need

Representing the local and global environments in atomic systems in a descriptive and efficient way has been an important research top in the chemistry, chemical engineering, and material science communities. Having good representations, or features, of chemical environments has proven to be vital for building reliable machine learning (ML) models. These models can accurately predict properties of atomic systems, and in limited cases have even been used for discovering or designing new chemicals and materials [@zuo_accelerating_2021; @collins_accelerated_2017]. So far, scientists and researchers have designed featurization schemes like the atom-centered symmetry function (ACSF) [@BehlerParrinello], the smooth overlap of atomic positions (SOAP) [@SOAP] and the Gaussian Momentum [@Gaussian_Momentum] schemes. More recently, graph representation and ML models based on them (e.g. MEGNet [@MEGNet], CGCNN [@CGCNN]) have been successful. Gaussian Multipole, or GMP [@GMP], is a recently developed scheme of featurizing local chemical environments, i.e. the chemical characteristics of spaces near individual atoms in molecules and crystal structures. GMP approximates underlying local electronic environments (e.g. approximated distribution of local electron cloud) using multipole expansion, the theory of which is explicated in a prior publication [@GMP]. The featurization scheme is flexible, depending only on prior assumptions of atomic identity and position, and it is therefore applicable to various atomic systems (molecules, nano-particles, periodic crystals, etc.) in which atomic arrangements are known. Feature computation is fast, and the representation accuracy is systematically improvable. Moreover, thanks to the deep connection between the GMP features and physics, we have previously shown that ML models based on these features are transferable [@GMP]. With these characteristics, GMP featurization could be useful to a broad audience for future research in chemistry and materials science. Therefore, having lightweight, reliable, and open software that can calculate these features in a fast and accurate way is desirable.

# Overview

The GMP-Featurizer package is mainly written in C++ and Python. C++ is used for the underlying computation module for speed, and Python used for an intuitive and readable API for scientists and researchers in the community to use. Two of the most widely used Python libraries in chemistry and materials science are ASE [@ase-paper] and pymatgen [@pymatgen]. The gmp-featurizer library provides APIs that can read atomic systems as input from these libraries, though they are not required as dependencies. This allows the users to leverage the I/O functionalities of these two packages to get data into the format that's accepted by GMP-Featurizer. On top of that, Ray is used to parallelize the feature computation, and the parallelization efficiency is close to 100\%.  Overall, this package is designed to be lightweight, easy to use, fast and accurate.

The main inputs of the workflow are a Python dictionary that contains the necessary hyperparameters for defining the desired GMP feature set, and a list of atomic systems that needs to be featured. The native way of defining atomic systems is simply a Python dictionary that contains information like `lattice vectors`, `atom positions`, `atom types`, etc. As mentioned, the package also supports both ASE Atoms and pymatgen Structure objects with pre-defined converters. This capability is extensible to other formats with custom-made converters. It also supports the featurization of disordered atomic structures, which is unsupported by many popular featurization methods. Please refer to the examples for more details. The output is simply a list of dictionaries containing the resulting features, and their derivatives if requested, for the atomic structures.

By default, the package computes GMP features at each atom position, but it can also be used to compute the features at any set of reference points inside the atomic system by providing a list of the positions of interest for each atomic structure. Users can also specify the number of cores for parallel computing. Noticeably, the GMP-featurizer package, unlike common local chemical environment featurization schemes (including the original implementation of GMP in AmpTorch), eliminates the need to define a hard distance cutoff for determining the local environment to be featurized. Instead, users simply need to state the desired accuracy of the features. GMP-featurizer then autonomously identifies neighboring atoms that will make meaningful contributions to each GMP feature relative to the set accuracy level, and computes the features. This change is significant, because GMP featurizes the underlying electron density of chemical systems, which can greatly vary between elements and hence influence the GMP features differently. For instance, a distant Mercury atom may have a larger impact on GMP features than a closer Hydrogen atom. Additionally, different GMP features respond differently to changes in atomic positions, making a universal cutoff inefficient and potentially leading to a lot of unnecessary computations (when the cutoff is set too high) or inaccurate features (when set too low). GMP-featurizer's new adaptive scheme circumvents these issues. It not only enhances efficiency and precision but also simplifies the user experience.

Moreover, computed results can be cached locally for convenient reprocessing of datasets, e.g. after augmentation or modification. Two sets of standard pseudopotential files are also provided, which are necessary to specify GMP feature sets, but may be difficult to collect from either commercial or open-source density functional theory systems. Lastly, a series of tutorials are provided in the repository to help users with quick starting and understanding the various features of the codebase.


# Acknowledgements

This work was supported by the Energy and Materials Division of the Toyota Research Institute. The authors acknowledge Jens Hummelshøj for helpful discussions regarding unit-testing frameworks.

# References