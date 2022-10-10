import numpy as np
import ray
from ray.util import ActorPool
from .base_feature import BaseFeature
from .GMP import GMP

import h5py
from tqdm import tqdm
from .util import get_hash, list_symbols_to_indices, validate_image


def calculate_GMP_features(
    feature_setup,
    elements,
    images,
    ref_positions_list=None,
    calc_derivatives=False,
    save_features=False,
    verbose=False,
    cores=1,
):

    # images_feature_list = []

    if ref_positions_list is None:
        ref_positions_list = [image.get_positions() for image in images]

    assert len(images) == len(ref_positions_list)

    # if save is true, create directories if not exist
    # feature._setup_feature_database(save_features=save_features)

    if cores <= 1:
        feature = GMP(feature_setup, elements)
        images_feature_list = []
        for image, ref_positions in tqdm(
            zip(images, ref_positions_list),
            total=len(images),
            desc="Computing features",
            disable=not verbose,
        ):
            temp_image_dict = feature._calculate_single_image(
                image,
                ref_positions,
                calc_derivatives,
                save_features,
            )
            images_feature_list.append(temp_image_dict)
        return images_feature_list

    elif cores > 1:

        remote_feature_actor = ray.remote(GMP)
        length = len(images)
        calc_deriv_list = [calc_derivatives] * length
        save_features_list = [save_features] * length
        args = zip(images, ref_positions_list, calc_deriv_list, save_features_list)

        ray.init(num_cpus=cores)
        actors = [
            remote_feature_actor.remote(feature_setup, elements) for _ in range(cores)
        ]
        pool = ActorPool(actors)
        poolmap = pool.map(
            lambda a, v: a._calculate_single_image.remote(v[0], v[1], v[2], v[3]), args
        )
        images_feature_list = [a for a in tqdm(poolmap, total=length)]

        ray.shutdown()
        return images_feature_list

        # length = len(images)
        # obj_ids = [remote_actor.remote(feature, image, ref_positions, calc_derivatives, save_features) for image, ref_positions in zip(images, ref_positions_list)]

        # for x in tqdm(to_iterator(obj_ids), total=length):
        #     pass

        # from .util import istarmap
        # import multiprocessing.pool as mpp
        # from multiprocessing import Pool

        # mpp.Pool.istarmap = istarmap

        # length = len(images)
        # calc_deriv_list = [calc_derivatives] * length
        # save_features_list = [save_features] * length
        # args = zip(images, ref_positions_list, calc_deriv_list, save_features_list)
        # images_feature_list = []
        # with Pool(cores) as p:
        #     for temp_image_dict in tqdm(
        #         p.istarmap(self._calculate_single_image, args), total=length
        #     ):
        #         images_feature_list.append(temp_image_dict)

        # return images_feature_list

    else:
        raise ValueError


class GMPFeaturizer:
    def __init__(
        self,
        GMPs,
        elements,
        calc_derivatives=False,
        save_features=False,
        verbose=True,
        # cores=1,
    ):

        # self.feature = GMP(GMPs, elements)
        self.feature_setup = GMPs
        self.elements = elements
        self.calc_derivatives = calc_derivatives
        self.save_features = save_features
        # self.cores = cores
        self.verbose = verbose

        # self.element_list = self.feature._get_element_list()
        self.features_ready = False

    def prepare_features(self, images, ref_positions_list=None, cores=1):

        if ref_positions_list == None:
            ref_positions_list = [image.get_positions() for image in images]

        # self.calculated_features_list = self.feature.calculate(
        self.calculated_features_list = calculate_GMP_features(
            self.feature_setup,
            self.elements,
            images,
            ref_positions_list,
            calc_derivatives=self.calc_derivatives,
            save_features=self.save_features,
            cores=cores,
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
