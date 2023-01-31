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

# setup featurizer/calcualtor 
sigmas = [0.2, 0.4, 0.6, 0.8, 1.0]
GMPs = {
    "GMPs": {   "orders": [-1, 0, 1, 2, 3], "sigmas":sigmas   },
    "psp_path": "/home/raylei/Documents/GMP-featurizer/pseudopotentials/NC-SR.gpsp",
    "square":False,
    "solid_harmonics": True,
}

featurizer = GMPFeaturizer(GMPs=GMPs)


# calculate features
result = featurizer.prepare_features(images)

# access data
features = [entry["features"] for entry in result]
feature_primes = [entry["feature_primes"] for entry in result]
```