# Copyright Toyota Research Institute 2023
"""
Module for defining the base class of features
"""

import os
from abc import ABC, abstractmethod
import h5py
import numpy as np

from .util import get_hash


class BaseFeature(ABC):
    """
    Base class for features
    """

    def __init__(self):
        """Base initialization"""
        super().__init__()
        self.feature_database = "cache/features/"

        self.feature_type = "default"
        self.feature_setup_hash = "default"

        self.database_initialized = False

    @abstractmethod
    def calculate_features(self, image, params_set, calculate_derivatives=True):
        """
        abstract method: calculating features
        """
        pass

    @abstractmethod
    def get_feature_setup_hash(self):
        """
        abstract method: getting hash for feature setup
        """
        pass

    @abstractmethod
    def save_feature_setup(self, filename):
        """
        abstract method: save feature setup to file
        """
        pass

    @abstractmethod
    def _prepare_feature_parameters(self):
        """
        abstract private method: prepare feature parameters based on config
        """
        pass

    def _calculate_single_image(
        self,
        image,
        ref_positions,
        calc_derivatives,
        calc_occ_derivatives,
        save_features,
        idx=0,
    ):
        """
        Control method for computing the features for a given image object

        Parameters
        ----------
        image : Dict
            Information of a image object in the format readable by this package
        ref_positions : numpy.ndarray
            Size n*3, Cartesian coordinates of reference points that need
            to be featurized
        calc_derivatives : bool
            Whether to calcualte feature derivatives w.r.t. atom positions
        calc_occ_derivatives : bool
            Whether to calcualte feature derivatives w.r.t. occupancies
        save_features : bool
            Whether to write resulting features to database
        idx : int (default 0)
            idx of the image, used to ordering

        Return
        ----------
        temp_image_dict : Dict
            Dictionary containing computed features, and derivatives if asked
        idx : int
            Needed for ordering 
        """

        # print("start")
        ref_positions = np.array(ref_positions)
        # validate_image(image, ref_positions)

        # if save, then read/write from db as needed
        if save_features:
            self._setup_feature_database(save_features=save_features)
            image_hash = get_hash(image, ref_positions)
            image_db_filename = "{}/{}.h5".format(
                self.desc_feature_database_dir, image_hash
            )
            try:
                temp_image_dict = self._compute_features(
                    image,
                    ref_positions,
                    image_db_filename,
                    calc_derivatives=calc_derivatives,
                    calc_occ_derivatives=calc_occ_derivatives,
                    save_features=save_features,
                )
            except Exception:
                print(
                    "File {} not loaded properly\nProceed to compute in run-time".format(
                        image_db_filename
                    )
                )
                temp_image_dict = self._compute_features_nodb(
                    image,
                    ref_positions,
                    calc_derivatives=calc_derivatives,
                    calc_occ_derivatives=calc_occ_derivatives,
                    save_features=save_features,
                )

        # if not save, compute fps on-the-fly
        else:
            temp_image_dict = self._compute_features_nodb(
                image,
                ref_positions,
                calc_derivatives=calc_derivatives,
                calc_occ_derivatives=calc_occ_derivatives,
                save_features=save_features,
            )

        return temp_image_dict, idx

    def _compute_features(
        self,
        image,
        ref_positions,
        image_db_filename,
        calc_derivatives,
        calc_occ_derivatives,
        save_features,
    ):
        """
        Computing the features for a given image object with database
        """
        with h5py.File(image_db_filename, "a") as db:
            image_dict = {}

            try:
                current_snapshot_grp = db[str(0)]
            except Exception:
                current_snapshot_grp = db.create_group(str(0))

            if calc_derivatives:
                try:
                    size_info = np.array(current_snapshot_grp["size_info"])
                    features = np.array(current_snapshot_grp["features"])
                    feature_primes_val = np.array(
                        current_snapshot_grp["feature_primes_val"]
                    )
                    feature_primes_row = np.array(
                        current_snapshot_grp["feature_primes_row"]
                    )
                    feature_primes_col = np.array(
                        current_snapshot_grp["feature_primes_col"]
                    )
                    feature_primes_size = np.array(
                        current_snapshot_grp["feature_primes_size"]
                    )
                    if calc_occ_derivatives:
                        feature_occ_primes = np.array(
                            current_snapshot_grp["feature_occ_primes"]
                        )
                except Exception:
                    (
                        size_info,
                        features,
                        feature_primes_val,
                        feature_primes_row,
                        feature_primes_col,
                        feature_primes_size,
                        feature_occ_primes,
                    ) = self.calculate_features(
                        image["cell"],
                        image["pbc"],
                        image["atom_positions"],
                        image["atom_symbols"],
                        image["occupancies"],
                        ref_positions,
                        calc_derivatives=calc_derivatives,
                        calc_occ_derivatives=calc_occ_derivatives,
                    )

                    if save_features:
                        current_snapshot_grp.create_dataset("size_info", data=size_info)
                        current_snapshot_grp.create_dataset("features", data=features)
                        current_snapshot_grp.create_dataset(
                            "feature_primes_val", data=feature_primes_val
                        )
                        current_snapshot_grp.create_dataset(
                            "feature_primes_row", data=feature_primes_row
                        )
                        current_snapshot_grp.create_dataset(
                            "feature_primes_col", data=feature_primes_col
                        )
                        current_snapshot_grp.create_dataset(
                            "feature_primes_size", data=feature_primes_size
                        )
                        if calc_occ_derivatives:
                            current_snapshot_grp.create_dataset(
                                "feature_occ_primes", data=feature_occ_primes
                            )

                image_dict["features"] = features
                # image_dict["num_features"] = size_info[2]
                feature_prime_dict = {}
                feature_prime_dict["row"] = feature_primes_row
                feature_prime_dict["col"] = feature_primes_col
                feature_prime_dict["val"] = feature_primes_val
                image_dict["feature_primes"] = feature_prime_dict
                if calc_occ_derivatives:
                    image_dict["feature_occ_primes"] = feature_occ_primes

            else:
                # features = np.array(current_snapshot_grp["features"])
                try:
                    # size_info = np.array(current_snapshot_grp["size_info"])
                    features = np.array(current_snapshot_grp["features"])
                    if calc_occ_derivatives:
                        feature_occ_primes = np.array(
                            current_snapshot_grp["feature_occ_primes"]
                        )
                except Exception:
                    (
                        size_info,
                        features,
                        _,
                        _,
                        _,
                        _,
                        feature_occ_primes,
                    ) = self.calculate_features(
                        image["cell"],
                        image["pbc"],
                        image["atom_positions"],
                        image["atom_symbols"],
                        image["occupancies"],
                        ref_positions,
                        calc_derivatives=calc_derivatives,
                        calc_occ_derivatives=calc_occ_derivatives,
                    )

                    if save_features:
                        current_snapshot_grp.create_dataset("size_info", data=size_info)
                        current_snapshot_grp.create_dataset("features", data=features)
                        if calc_occ_derivatives:
                            current_snapshot_grp.create_dataset(
                                "feature_occ_primes", data=feature_occ_primes
                            )

                image_dict["features"] = features
                if calc_occ_derivatives:
                    image_dict["feature_occ_primes"] = feature_occ_primes

        return image_dict

    def _compute_features_nodb(
        self,
        image,
        ref_positions,
        calc_derivatives,
        calc_occ_derivatives,
        save_features,
    ):
        """
        Computing the features for a given image object without database
        """
        image_dict = {}

        if calc_derivatives:

            (
                size_info,
                features,
                feature_primes_val,
                feature_primes_row,
                feature_primes_col,
                feature_primes_size,
                feature_occ_primes,
            ) = self.calculate_features(
                image["cell"],
                image["pbc"],
                image["atom_positions"],
                image["atom_symbols"],
                image["occupancies"],
                ref_positions,
                calc_derivatives=calc_derivatives,
                calc_occ_derivatives=calc_occ_derivatives,
            )

            image_dict["features"] = features
            # image_dict["num_features"] = size_info[2]
            feature_prime_dict = {}
            feature_prime_dict["size"] = feature_primes_size
            feature_prime_dict["row"] = feature_primes_row
            feature_prime_dict["col"] = feature_primes_col
            feature_prime_dict["val"] = feature_primes_val
            image_dict["feature_primes"] = feature_prime_dict
            if calc_occ_derivatives:
                image_dict["feature_occ_primes"] = feature_occ_primes

        else:
            (
                size_info,
                features,
                _,
                _,
                _,
                _,
                feature_occ_primes,
            ) = self.calculate_features(
                image["cell"],
                image["pbc"],
                image["atom_positions"],
                image["atom_symbols"],
                image["occupancies"],
                ref_positions,
                calc_derivatives=calc_derivatives,
                calc_occ_derivatives=calc_occ_derivatives,
            )

            image_dict["features"] = features
            if calc_occ_derivatives:
                image_dict["feature_occ_primes"] = feature_occ_primes

        return image_dict

    # def _feature_prime_element_row_index_to_image_row_index(
    #     self, original_rows, index_arr, num_desc, num_desc_max
    # ):
    #     atom_indices_for_specific_element, desc_indices = np.divmod(
    #         original_rows, num_desc
    #     )

    #     atom_indices_in_image = index_arr[atom_indices_for_specific_element]

    #     new_row = atom_indices_in_image * num_desc_max + desc_indices
    #     return new_row

    def _setup_feature_database(self, save_features):
        """
        Initialize database if involved
        """
        if save_features:
            if self.database_initialized:
                return
            self.get_feature_setup_hash()
            self.desc_type_database_dir = os.path.join(
                self.feature_database, self.feature_type
            )

            self.desc_feature_database_dir = os.path.join(
                self.desc_type_database_dir, self.feature_setup_hash
            )

            os.makedirs(self.feature_database, exist_ok=True)
            os.makedirs(self.desc_type_database_dir, exist_ok=True)
            os.makedirs(self.desc_feature_database_dir, exist_ok=True)
            self.database_initialized = True
            feature_setup_filename = "feature_setup.dat"
            feature_setup_path = "{}/{}".format(
                self.desc_feature_database_dir, feature_setup_filename
            )
            self.save_feature_setup(feature_setup_path)

    # def get_element_list(self):
    #     """
    #     return the list of elements
    #     """
    #     return self.elements
