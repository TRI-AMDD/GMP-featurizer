# GMP-featurizer
![Testing - main](https://github.com/TRI-AMDD/GMP-featurizer/workflows/Testing%20-%20main/badge.svg)
![Linting](https://github.com/TRI-AMDD/GMP-featurizer/workflows/Linting/badge.svg)

This package is used to efficiently and accurately compute the GMP features and their derivatives for any chemical systems. The computation is also parallelized via Ray.

The details of the theory behind the Gaussian Multipole descriptors can be found in the [original paper](https://pubs.acs.org/doi/10.1021/acs.jpclett.2c02100)
or in its [arxiv version](https://arxiv.org/abs/2102.02390)

Part of the code of this package is based on the [AmpTorch package](https://github.com/ulissigroup/amptorch)

## Installation
To install this package, simply clone this repo, 
```
git clone https://github.com/TRI-AMDD/GMP-featurizer
cd GMP-featurizer
```

Then install the requirements and the package itself
```
pip install -r requirements.txt
pip install -e .
```

## Basic usage

#### Please refer to the example notebooks for better and detailed tutorials

#### Import modules and load data
```
import numpy as np
from GMPFeaturizer import GMPFeaturizer, ASEAtomsConverter, PymatgenStructureConverter
from ase.io import read as aseread

# Loading cif file as a ase atoms object
image = aseread("./test.cif") 
# The input to the featurizer should be a non-empty list
images = [image]

# initialize the converter, in this case it's the converter for ASE atoms objects
# There is also a pre-existing converter for pymatgen Structure objects as well
converter = ASEAtomsConverter()
# converter = PymatgenStructureConverter()
```

#### Setup the featurizer
The list of features is the Cartesian product of orders and sigams (except for order -1, which correspond just local electron density, so different simgas does not matter. Thus, there is only one feature for order -1). 

With this setting, the list of features are

[(-1, 0), (0, 0.1), (0, 0.2), (0, 0.3), (1, 0.1), (1, 0.2), (1, 0.3), (2, 0.1), (2, 0.2), (2, 0.3)]

where the first number is the order of the MCSH angular probe, and the second number is the sigma of the Gaussian radial probe 
```
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
    #no need to change if you are not considering the feature derivatives
    # "square": False, 
}

featurizer = GMPFeaturizer(GMPs=GMPs, calc_derivatives=True, verbose=True)
```
Set calc_derivatives=True if you want to get the feature derivatives w.r.t. atom positions, which are stored in the form of sparse matrices.



#### Calculate features and access data
Use the "cores" argument to change the number of cores for parallelization. Also converted needed to be specified,
```
result = featurizer.prepare_features(images, cores=5, converter=converter)

features = [entry["features"] for entry in result]
feature_primes = [entry["feature_primes"] for entry in result]
```

#### Specifying the list of GMP features
It's also possible to manually specify the list of GMP features to be computed, instead of specifying orders and sigmas.
```
GMPs = {
    "GMPs_detailed_list": [(-1,0), (0, 0.1), (0, 0.2), (0, 0.3), (1, 0.2), (1, 0.3), (2, 0.3)],
    "psp_path": "./NC-SR.gpsp", # path to the pseudo potential file
    "overlap_threshold": 1e-16, # basically the accuracy of the resulting features
    # "square": False, # whether the features are squared, no need to change if you are not get the feature derivatives
}
```

#### Whole Script
```
import numpy as np
from GMPFeaturizer import GMPFeaturizer, ASEAtomsConverter, PymatgenStructureConverter

# load data
from ase.io import read as aseread
image = aseread("./test.cif") 
images = [image]

converter = ASEAtomsConverter()
# converter = PymatgenStructureConverter()

# setup featurizer
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
    #no need to change if you are not considering the feature derivatives
    # "square": False, 
}
featurizer = GMPFeaturizer(GMPs=GMPs, calc_derivatives=True, verbose=True)



# calculate features
result = featurizer.prepare_features(images, cores=5, converter=converter)

# access data
features = [entry["features"] for entry in result]
feature_primes = [entry["feature_primes"] for entry in result]
```

#### Save calculated feature to / load calculated feature from local folder
Simply set "save_features=True" when calling the prepare_features function.

The path to the local database is set when initializing the featurizer
```
featurizer = GMPFeaturizer(GMPs=GMPs, calc_derivatives=False, feature_database="cache/features/")
features = featurizer.prepare_features(images, cores=5, save_features=True, converter=converter)
```

#### License
Apache 2.0

Copyright 2023 Toyota Research Institute