# GMP-featurizer
This package is used to efficiently and accurately compute the GMP features and their derivatives for any chemical systems. The computation is also parallelized via Ray.

## Basic usage

#### Import modules and load data
```
import numpy as np
from GMPFeaturizer.GMP import GMP
from ase.io import read as aseread

image = aseread("./test.cif")
images = [image] # it should be a non-empty list
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
    "psp_path": "./NC-SR.gpsp", # path to the pseudo potential file
    "overlap_threshold": 1e-16, # basically the accuracy of the resulting features
    # "square": False, # whether the features are squared, no need to change if you are not get the feature derivatives
}

featurizer = GMPFeaturizer(GMPs=GMPs, calc_derivatives=True)
```
set calc_derivatives=True if you want to get the feature derivatives w.r.t. atom positions, which are stored in the form of sparse matrices


#### Calculate features and access data
Use the "cores" argument to change the number of cores for parallelization
```
result = featurizer.prepare_features(images, cores=5)

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
from GMPFeaturizer.GMP import GMP

# load data
from ase.io import read as aseread
image = aseread("./test.cif")
images = [image] 

# setup featurizer
GMPs = {
    "GMPs": {   
        "orders": [-1, 0, 1, 2], 
        "sigmas": [0.1, 0.2, 0.3]   
    },
    "psp_path": "./NC-SR.gpsp", # path to the pseudo potential file
    "overlap_threshold": 1e-16, # basically the accuracy of the resulting features
    # "square": False, # whether the features are squared, no need to change if you are not get the feature derivatives
}
featurizer = GMPFeaturizer(GMPs=GMPs, calc_derivatives=True)



# calculate features
result = featurizer.prepare_features(images, cores=5)

# access data
features = [entry["features"] for entry in result]
feature_primes = [entry["feature_primes"] for entry in result]
```