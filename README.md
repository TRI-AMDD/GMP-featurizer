# GMP-featurizer
Feature calculator for GMP features

## basic usage
```
import numpy as np
from GMPFeaturizer.GMPOrderNorm import GMPOrderNorm

# load data
from ase.io import read as aseread
image = aseread("./test.cif")
images = [image]

#GMP hyperparameters
GMP_order = 3
nsigmas = 5
width = 1.0
elements = ["Ce", "Zr", "Ti", "Mg", "O"]


# setup featurizer/calcualtor 
sigmas = np.linspace(0.0,width,nsigmas+1,endpoint=True)[1:]
GMPs = {
            "GMPs": {   "orders": [-1]+list(range(GMP_order+1)), "sigmas":sigmas   },
            "atom_gaussians": {x:"<path_to_psp>/{}_pseudodensity.g".format(x) for x in elements},
            #"cutoff": 5.0,
            "square":False,
            "solid_harmonics": True,
}

featurizer = GMPOrderNorm(GMPs=GMPs, elements=elements)


# calculate features
result = featurizer.calculate(images, calc_derivatives=True)

# access data
features = [entry["features"] for entry in result]
feature_primes = [entry["feature_primes"] for entry in result]
```