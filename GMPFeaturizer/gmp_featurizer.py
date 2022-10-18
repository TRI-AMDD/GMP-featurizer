import os
import numpy as np
import ase
import ray
from ray.util import ActorPool
from .base_feature import BaseFeature
from .GMP import GMP

import h5py
from tqdm import tqdm
from .util import get_hash, list_symbols_to_indices  # , validate_image
from .converters import ASEAtomsConverter, PymatgenStructureConverter


class GMPFeaturizer:
    def __init__(
        self, GMPs, elements, calc_derivatives=False, verbose=True,
    ):

        # self.feature = GMP(GMPs, elements)
        self.feature_setup = GMPs
        self.elements = elements
        self.calc_derivatives = calc_derivatives
        # self.save_features = save_features
        # self.cores = cores
        self.verbose = verbose

    def prepare_features(
        self,
        image_objects,
        ref_positions_list=None,
        cores=1,
        save_features=False,
        converter=None,
    ):

        images = self._convert_validate_image_objects(image_objects, converter)

        if ref_positions_list == None:
            ref_positions_list = [image["atom_positions"] for image in images]

        assert len(images) == len(ref_positions_list)

        if cores == 0:
            cores = os.cpu_count()

        # self.calculated_features_list = self.feature.calculate(
        # calculated_features_list = calculate_GMP_features(
        calculated_features_list = self._calculate_GMP_features(
            images,
            ref_positions_list=ref_positions_list,
            calc_derivatives=self.calc_derivatives,
            save_features=save_features,
            cores=cores,
            verbose=self.verbose,
        )

        return calculated_features_list

    def _convert_validate_image_objects(self, image_objects, converter):
        if converter is None:
            if isinstance(image_objects[0], ase.Atoms):
                converter = ASEAtomsConverter()
                images = converter.convert(image_objects)
            else:
                images = image_objects
        else:
            images = converter(image_objects)

        for image in images:
            assert isinstance(image, dict)
            assert "pbc" in image
            assert "atom_positions" in image
            assert "atom_symbols" in image
            assert "cell" in image
            if "occupancies" not in image:
                image["occupancies"] = np.array(
                    [1.0 for _ in range(len(image["atom_symbols"]))]
                )
            assert len(image["atom_positions"]) == len(image["atom_symbols"])
            assert len(image["atom_positions"]) == len(image["occupancies"])

        return images

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

    def _calculate_GMP_features(
        self,
        images,
        ref_positions_list=None,
        calc_derivatives=False,
        save_features=False,
        verbose=False,
        cores=1,
    ):
        # if ref_positions_list is None:
        #     ref_positions_list = [image.get_positions() for image in images]

        # assert len(images) == len(ref_positions_list)

        # if save is true, create directories if not exist
        # feature._setup_feature_database(save_features=save_features)

        if cores <= 1:
            feature = GMP(self.feature_setup, self.elements)
            images_feature_list = []
            for image, ref_positions in tqdm(
                zip(images, ref_positions_list),
                total=len(images),
                desc="Computing features",
                disable=not verbose,
            ):
                temp_image_dict = feature._calculate_single_image(
                    image, ref_positions, calc_derivatives, save_features,
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
                remote_feature_actor.remote(self.feature_setup, self.elements)
                for _ in range(cores)
            ]
            pool = ActorPool(actors)
            poolmap = pool.map(
                lambda a, v: a._calculate_single_image.remote(v[0], v[1], v[2], v[3]),
                args,
            )
            images_feature_list = [a for a in tqdm(poolmap, total=length)]

            ray.shutdown()
            return images_feature_list
