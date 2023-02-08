# GMP-featurizer
Feature calculator for GMP features

## basic usage
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

# The list of features is the Cartesian product of orders and sigams (except for order -1, which correspond just local electron density, 
# so different simgas does not matter. Thus, there is only one feature for order -1). 
# With this setting, the list of features are
# [(-1, 0), (0, 0.1), (0, 0.2), (0, 0.3), (1, 0.1), (1, 0.2), (1, 0.3), (2, 0.1), (2, 0.2), (2, 0.3)]
#where the first number is the order of the MCSH angular probe, and the second number is the sigma of the Gaussian radial probe 


# set calc_derivatives=True if you want to get the feature derivatives w.r.t. atom positions,
# which are stored in the form of sparse matrices
featurizer = GMPFeaturizer(GMPs=GMPs, calc_derivatives=True)



# calculate features
result = featurizer.prepare_features(images)

# access data
features = [entry["features"] for entry in result]
feature_primes = [entry["feature_primes"] for entry in result]
```

You can also manually specify the features to be computed:
```
GMPs = {
    "GMPs_detailed_list": [(-1,0), (0, 0.1), (0, 0.2), (0, 0.3), (1, 0.2), (1, 0.3), (2, 0.3)],
    "psp_path": "./NC-SR.gpsp", # path to the pseudo potential file
    "overlap_threshold": 1e-16, # basically the accuracy of the resulting features
    # "square": False, # whether the features are squared, no need to change if you are not get the feature derivatives
}
```