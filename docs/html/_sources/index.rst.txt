.. GMP-featurizer documentation master file, created by
   sphinx-quickstart on Thu Jun  8 17:55:40 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. toctree::
   :maxdepth: 2


==============
GMP-featurizer
==============

Introduction
==========================================

This package is used to efficiently compute the GMP features and their derivatives for any chemical system. The computation is also parallelized via Ray.  The details of the theory behind the Gaussian Multipole descriptors can be found in the `original paper <https://pubs.acs.org/doi/10.1021/acs.jpclett.2c02100>`_ or in its `arxiv version <https://arxiv.org/abs/2102.02390>`_. Part of the code of this package is based on the `AmpTorch package <https://github.com/ulissigroup/amptorch>`_, which is gratefully acknowledged.

For a quickstart tutorial - we recommend reading the :ref:`Basic usage` section below. For more detailed examples, please browser the `examples folder in the github repository <https://github.com/TRI-AMDD/GMP-featurizer/tree/docs/examples>`_.

To browse the API, associated functions, please refer to the :ref:`Module index.<modindex>`.

Installation
==============
To install this package, clone the `repo <https://github.com/TRI-AMDD/GMP-featurizer>`_ using git.

.. code-block:: console

   git clone https://github.com/TRI-AMDD/GMP-featurizer
   cd GMP-featurizer

Then install the requirements and the package itself

.. code-block:: console

   pip install -r requirements.txt
   pip install -e .

Basic usage
=============

Shown below is a basic tutorial. Please refer to the `example notebooks <https://github.com/TRI-AMDD/GMP-featurizer/tree/docs/examples>`_ for detailed tutorials. An example "cif" file is provided in the "examples" directory

.. code-block:: python

   # Import modules and load data
   import numpy as np
   from GMPFeaturizer import GMPFeaturizer, ASEAtomsConverter, PymatgenStructureConverter
   from ase.io import read as aseread

   # Loading cif file as a ase atoms object
   image = aseread("./examples/test.cif")

   # The input to the featurizer should be a non-empty list
   images = [image]

   # initialize the converter, in this case it's the converter for ASE atoms objects
   # There is also a pre-existing converter for pymatgen Structure objects as well

   converter = ASEAtomsConverter()
   converter = PymatgenStructureConverter()

   ### Setup the featurizer
   # The list of features is the Cartesian product of orders and sigams
   # (except for order -1, which correspond just local electron density,
   # so different sigmas does not matter. Thus, there is only one feature for order -1).
   # With this setting, the list of features are
   # [(-1, 0), (0, 0.1), (0, 0.2), (0, 0.3), (1, 0.1),
   #  (1, 0.2), (1, 0.3), (2, 0.1), (2, 0.2), (2, 0.3)]
   # where the first number is the order of the MCSH angular probe,
   # and the second number is the sigma of the Gaussian radial probe
   GMPs = {
       "GMPs": {
           "orders": [-1, 0, 1, 2],
           "sigmas": [0.1, 0.2, 0.3]
       },
       # path to the pseudo potential file
       "psp_path": "<path>/NC-SR.gpsp",

       # basically the accuracy of the resulting features
       "overlap_threshold": 1e-16,

       # whether the features are squared,
       # no need to change if you are not considering the feature derivatives
       # "square": False,
   }

   featurizer = GMPFeaturizer(GMPs=GMPs, converter=converter, calc_derivatives=True, verbose=True)

   # Set calc_derivatives=True if you want to get the feature derivatives w.r.t. atom positions, which are stored in the form of sparse matrices.

   ### Calculate features and access data
   # Use the "cores" argument to change the number of cores for parallelization. Also converted needed to be specified,

   result = featurizer.prepare_features(images, cores=5)
   features = [entry["features"] for entry in result]
   feature*primes = [entry["feature*primes"] for entry in result]


Specifying the list of GMP features
###################################

It's also possible to manually specify the list of GMP features to be computed, instead of specifying orders and sigmas.

.. code-block:: python

   GMPs = {
       "GMPs_detailed_list": [(-1,0), (0, 0.1), (0, 0.2), (0, 0.3), (1, 0.2), (1, 0.3), (2, 0.3)],
       "psp_path": "./NC-SR.gpsp", # path to the pseudo potential file
       "overlap_threshold": 1e-16, # basically the accuracy of the resulting features
       # "square": False, # whether the features are squared, no need to change if you are not get the feature derivatives
   }

Copyright 2023 Toyota Research Institute



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
