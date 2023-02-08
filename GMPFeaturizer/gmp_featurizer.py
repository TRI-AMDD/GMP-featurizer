import os
import numpy as np
import ase
import ray
from ray.util import ActorPool
from .GMP import GMP

from tqdm import tqdm
from .converters import ASEAtomsConverter


class GMPFeaturizer:
    """ object for computing GMP features of checmical systems """
    def __init__(
        self,
        GMPs,
        feature_database="cache/features/",
        calc_derivatives=False,
        calc_occ_derivatives=False,
        verbose=True,
    ):
        """
        Parameters
        ----------
        GMPs : dict
            configuration dictionary for the GMP feature set
        feature_database : str, optional (default: "cache/features/")
            path to the database for storing calculated features
        calc_derivative : bool (default: False)
            Whether to calculate feature derivatives w.r.t. atom positions
        calc_occ_derivative : bool (default: False)
            Whether to calculate feature derivatives w.r.t. atom occupancies
        verbose : bool (default: True)
            Whether to print information
        """

        self.feature_setup = GMPs
        self.feature_database = feature_database
        self.calc_derivatives = calc_derivatives
        self.calc_occ_derivatives = calc_occ_derivatives
        self.verbose = verbose

    def prepare_features(
        self,
        image_objects,
        ref_positions_list=None,
        cores=1,
        save_features=False,
        converter=None,
    ):
        """
        computing features with given list of image objects

        Parameters
        ----------
        image_objects : list
            list of image objects
        ref_positions_list : list
            list of lists that contains for positions for computing
            GMP features
        cores : int (default: 1)
            number of cores for computing the features
        save_features : bool (default: False)
            whether to save features to database
        converter : converter_object (default: None)
            converter to convert image objects to objects that 
            can be read. Examples can be found and imported from 
            GMPFeaturizer.ASEAtomsConverter
            GMPFeaturizer.PymatgenStructureConverter

        Return
        ----------
        calculated_features_list : list
            list of dicts that store the computed features and derivatives 
        """

        images = self._convert_validate_image_objects(image_objects, converter)

        if ref_positions_list is None:
            ref_positions_list = [image["atom_positions"] for image in images]

        assert len(images) == len(ref_positions_list)

        if cores == 0:
            cores = os.cpu_count()

        calculated_features_list = self._calculate_GMP_features(
            images,
            ref_positions_list=ref_positions_list,
            calc_derivatives=self.calc_derivatives,
            calc_occ_derivatives=self.calc_occ_derivatives,
            save_features=save_features,
            cores=cores,
            verbose=self.verbose,
        )

        return calculated_features_list

    def _convert_validate_image_objects(self, image_objects, converter):
        """
        Private method for validating the image objects
        """
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
        """
        PCA transform of features, not implemented yet
        """
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
        """
        Scale features, not implemented yet
        """
        raise NotImplementedError

    # TODO
    def get_stats(self):
        """
        Get statistical information of the features, not implemented yet
        """
        raise NotImplementedError

    def _calculate_GMP_features(
        self,
        images,
        ref_positions_list=None,
        calc_derivatives=False,
        calc_occ_derivatives=False,
        save_features=False,
        verbose=False,
        cores=1,
    ):
        """
        Private method for parallelized computation of the features
        """
        # if ref_positions_list is None:
        #     ref_positions_list = [image.get_positions() for image in images]

        # assert len(images) == len(ref_positions_list)

        # if save is true, create directories if not exist
        # feature._setup_feature_database(save_features=save_features)

        if cores <= 1:
            feature = GMP(self.feature_setup, self.feature_database)
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
                    calc_occ_derivatives,
                    save_features,
                )
                images_feature_list.append(temp_image_dict)
            return images_feature_list

        elif cores > 1:

            remote_feature_actor = ray.remote(GMP)
            length = len(images)
            calc_deriv_list = [calc_derivatives] * length
            calc_occ_deriv_list = [calc_occ_derivatives] * length
            save_features_list = [save_features] * length
            args = zip(
                images,
                ref_positions_list,
                calc_deriv_list,
                calc_occ_deriv_list,
                save_features_list,
            )

            ray.init(num_cpus=cores)
            actors = [
                remote_feature_actor.remote(self.feature_setup, self.feature_database)
                for _ in range(cores)
            ]
            pool = ActorPool(actors)
            poolmap = pool.map(
                lambda a, v: a._calculate_single_image.remote(
                    v[0], v[1], v[2], v[3], v[4]
                ),
                args,
            )
            images_feature_list = [a for a in tqdm(poolmap, total=length, disable=not verbose)]

            ray.shutdown()
            return images_feature_list
