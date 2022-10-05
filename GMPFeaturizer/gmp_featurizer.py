import numpy as np

from .base_feature import BaseFeature
from .GMPOrderNorm import GMPOrderNorm


class GMPFeaturizer:
    def __init__(
        self,
        GMPs,
        elements,
        calc_derivatives=False,
        save_features=False,
        verbose=True,
        cores=1,
    ):

        self.feature = GMPOrderNorm(GMPs, elements)
        self.calc_derivatives = calc_derivatives
        self.save_features = save_features
        self.cores = cores
        self.verbose = verbose

        self.element_list = self.feature._get_element_list()
        self.features_ready = False

    def prepare_features(self, images, ref_positions_list=None):

        if ref_positions_list == None:
            ref_positions_list = [image.get_positions() for image in images]

        self.calculated_features_list = self.feature.calculate(
            images,
            ref_positions_list,
            calc_derivatives=self.calc_derivatives,
            save_features=self.save_features,
            cores=self.cores,
            verbose=self.verbose,
        )

        self.features_ready = True

        return self.calculated_features_list

    def get_features(self):
        if not self.features_ready:
            print(
                "ERROR, features not calculated yet, please call prepare_features() function first"
            )
            return None

        else:
            return self.calculated_features_list

    # TODO
    def calculate_PCA(self, save_models=True, n_components=10, apply_PCA=True):
        raise NotImplementedError

    # TODO
    def calculate_scaling(
        self,
        separate_atomtypes=True,
        save_models=True,
        scaler_min=-1,
        scaler_max=1,
        apply_scaling=True,
    ):
        raise NotImplementedError

    # TODO
    def get_stats(self):
        raise NotImplementedError
