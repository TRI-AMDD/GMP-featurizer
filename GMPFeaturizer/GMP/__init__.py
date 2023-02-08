import hashlib
import math

import numpy as np
from scipy import sparse

from ..base_feature import BaseFeature
from ..constants import ATOM_SYMBOL_TO_INDEX_DICT
from ..util import (
    _gen_2Darray_for_ffi,
    list_symbols_to_indices,
    _gen_2Darray_for_ffi2,
    get_scaled_position,
)
from ._libgmp import ffi, lib


class GMP(BaseFeature):
    """
    Class for GMP feature
    """

    def __init__(
        self,
        GMPs,
        feature_database="cache/features/",
    ):
        """
        Parameters
        ----------
        GMPs: dict
            dictionary for the configuration of the features
            here are the parameters involved

            "GMPs" : dict
                dict for specifying which GMP features to be
                computed. it should have the following items:
                "orders" : list
                    list of orders of the intended GMP probes
                    note that order -1 means getting the local
                    electronic density
                "sigmas" : list
                    list of sigmas of the intended GMP probes
                    the final list of GMP probes(features) is the
                Cartesian product of the "orders" and "sigmas" lists
                except for order -1, which just correspond to local
                electronic density
            "GMPs_detailed_list" : list (default: None)
                list of GMP probes (features), e.g.
                [(-1, 0), (0, 0.1), (0, 0.2), ... (2, 0.5), (3, 0.3)...]
                Note that if this list is provided, "orders" and "sigmas"
                will be ignored.
            "psp_path" : str
                path to the GMP pseudopotential file (.gpsp file)
            "solid_harmonics" : bool (default: True)
                Whether to use solid harmonics
            "square" : book (default: False)
                Whether the resulting features are squared
            "custom_cutoff" : int (default: 4)
                How to calcualte the cutoff distance for feature computation
                different setting could mean very different computation speed
                of the features. No need to change if not farmiliar with this
                here are the possible settings
                -1 : manually specify cutoff distance
                0, 1, 2, 3, 4 : set cutoff based on the expected accuracy of the
                                resulting features, which is based on the overlap
                                between the probe functions and the gaussians in
                                the pseudo potentials of elements. These settings
                                are here mainly for accelerating the computation
                                of the GMP features without losing accuracy.
            "cutoff" : float
                manual setting of the cutof distance,
                needed when "custom_cutoff" is -1
            "overlap_threshold" : float (default: 1e-11)
                overlap threshold between the GMP probe function
                and elementalpseudopotential gaussians.
                Basically is the same as the expected accuracy of
                the resulting features. Needed when "custom_cutoff" is not -1

        feature_database : str (default: "cache/features/")
            path to the database of calculated features

        """
        super().__init__()
        self.feature_type = "GMP"
        self.feature_database = feature_database
        self.GMPs = GMPs
        self.custom_cutoff = self.GMPs.get("custom_cutoff", 4)

        self._load_psp_files()
        self._get_desc_list()
        self.overlap_threshold = self.GMPs.get("overlap_threshold", 1e-11)
        self._get_cutoffs()
        self._get_feature_setup()
        self._prepare_feature_parameters()

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, BaseFeature):
            if self.feature_type != other.feature_type:
                return False
            if self.feature_setup_hash != other.feature_setup_hash:
                return False
            return True
        return NotImplemented

        return cutoff

    # *******************************************************
    # related to cutoff calculation
    # *******************************************************

    def _get_C1_C2(self, A, B, alpha, beta):
        """
        private method for calculating C1 and C2
        """
        if alpha == 0:
            C1 = B
            C2 = -1.0 * beta
        else:
            temp = math.sqrt(math.pi / (alpha + beta))
            C1 = A * B * temp * temp * temp
            C2 = -1.0 * (alpha * beta / (alpha + beta))
        return C1, C2

    def _get_default_cutoff_element_gaussians(
        self, sigma, psp, max_num_gaussians, threshold=1e-8
    ):
        """
        private method for getting cutoffs based on
        the probe sigma and widest element gaussian
        """
        if sigma == 0:
            A = 0
            alpha = 0
        else:
            A = 1.0 / (sigma * math.sqrt(2.0 * math.pi))
            alpha = 1.0 / (2.0 * sigma * sigma)
        log_threshold = math.log(threshold)

        total_log_C1 = 0
        total_C2 = 0
        distances_list = np.zeros(max_num_gaussians)
        for i, gaussian in enumerate(psp):
            B, beta = gaussian
            C1, C2 = self._get_C1_C2(A, B, alpha, beta)
            C1 = abs(C1)
            temp_distance = math.sqrt((log_threshold - C1) / C2)
            distances_list[i] = temp_distance

        return distances_list.tolist()

    def _get_default_cutoff_single_element(self, sigma, psp, threshold=1e-8):
        """
        private method for getting cutoffs based on
        the probe sigma and widest element gaussian
        """
        if sigma == 0:
            A = 0
            alpha = 0
        else:
            A = 1.0 / (sigma * math.sqrt(2.0 * math.pi))
            alpha = 1.0 / (2.0 * sigma * sigma)
        log_threshold = math.log(threshold)

        total_log_C1 = 0
        total_C2 = 0
        distances_list = []
        for gaussian in psp:
            B, beta = gaussian
            C1, C2 = self._get_C1_C2(A, B, alpha, beta)
            C1 = abs(C1)
            temp_distance = math.sqrt((log_threshold - C1) / C2)
            distances_list.append(temp_distance)

        return np.max(distances_list).item()

    def _get_cutoffs(self):
        """
        get cutoffs based on setting
        """
        if self.custom_cutoff == -1:
            assert "cutoff" in self.GMPs
            self.cutoff = self.GMPs["cutoff"]
        elif self.custom_cutoff == 0:
            self._get_sigma_cutoffs()
        elif self.custom_cutoff == 1:
            self._get_sigma_cutoffs()
        elif self.custom_cutoff == 2:
            self._get_sigma_elemental_cutoffs()
        elif self.custom_cutoff == 3:
            self._get_order_sigma_elemental_cutoffs()
        elif self.custom_cutoff == 4:
            self._get_order_sigma_elemental_gaussian_cutoffs()

    def _get_order_sigma_elemental_cutoffs(self):
        """
        private method for getting cutoff distance based on probe order,
        probe gaussian sigma, and widest elemental gaussian
        """
        factors = {
            -1: 1,
            0: 1,
            1: 1,
            2: 1,
            3: 1e-1,
            4: 1e-2,
            5: 1e-4,
            6: 1e-5,
            7: 1e-7,
            8: 1e-9,
            9: 1e-11,
        }
        elemental_order_sigma_cutoffs = []
        order_sigma_index_dict = {}
        result = {}
        order_sigma_index = 0
        for order, sigma in self.desc_list:
            if (order, sigma) in order_sigma_index_dict:
                continue

            if order == -1:
                element_distances = []
                for element, psp in self.atomic_psp.items():
                    element_distance = self._get_default_cutoff_single_element(
                        0, psp, threshold=self.overlap_threshold
                    )
                    element_distances.append(element_distance)
                elemental_order_sigma_cutoffs.append(element_distances)
                order_sigma_index_dict[0] = order_sigma_index
                order_sigma_index += 1

            else:
                element_distances = []
                for element, psp in self.atomic_psp.items():
                    element_distance = self._get_default_cutoff_single_element(
                        sigma,
                        psp,
                        threshold=self.overlap_threshold * factors[order],
                    )
                    element_distances.append(element_distance)
                    # print(order, sigma, element, element_distance)
                elemental_order_sigma_cutoffs.append(element_distances)
                order_sigma_index_dict[(order, sigma)] = order_sigma_index
                order_sigma_index += 1

        self.order_sigma_index_dict = order_sigma_index_dict
        self.nsigmas = len(elemental_order_sigma_cutoffs)

        elemental_order_sigma_cutoffs = np.asarray(
            elemental_order_sigma_cutoffs, dtype=np.float64, order="C"
        )

        self.params_set["elemental_order_sigma_cutoffs"] = elemental_order_sigma_cutoffs
        # [[sigma1_element1_cutoff, sigma1_element2_cutoff,...]
        #  [sigma2_element1_cutoff, sigma2_element2_cutoff,...]
        #  ...
        #  [sigman_element1_cutoff, sigman_element2_cutoff,...]]
        self.params_set["elemental_order_sigma_cutoffs_p"] = _gen_2Darray_for_ffi(
            elemental_order_sigma_cutoffs, ffi
        )
        return

    def _get_order_sigma_elemental_gaussian_cutoffs(self):
        """
        private method for getting cutoff distance based on probe order,
        probe gaussian sigma, and each elemental gaussian
        """
        factors = {
            -1: 1,
            0: 1,
            1: 1,
            2: 1,
            3: 1e-1,
            4: 1e-2,
            5: 1e-4,
            6: 1e-5,
            7: 1e-7,
            8: 1e-9,
            9: 1e-11,
        }
        elemental_order_sigma_cutoffs = []
        elemental_order_sigma_gaussian_cutoffs = []

        self.max_num_gaussians = 0
        for psp in self.atomic_psp.values():
            if len(psp) > self.max_num_gaussians:
                self.max_num_gaussians = len(psp)

        order_sigma_index_dict = {}
        result = {}
        order_sigma_index = 0
        for order, sigma in self.desc_list:
            if (order, sigma) in order_sigma_index_dict:
                continue

            if order == -1:
                element_gaussian_distances = []
                element_distances = []
                for element, psp in self.atomic_psp.items():
                    gaussian_distances = self._get_default_cutoff_element_gaussians(
                        0, psp, self.max_num_gaussians, threshold=self.overlap_threshold
                    )
                    element_distance = np.max(gaussian_distances).item()
                    element_distances.append(element_distance)
                    element_gaussian_distances += gaussian_distances
                    # print(order, 0, element, element_distances)
                elemental_order_sigma_gaussian_cutoffs.append(
                    element_gaussian_distances
                )
                elemental_order_sigma_cutoffs.append(element_distances)
                order_sigma_index_dict[0] = order_sigma_index
                order_sigma_index += 1

            else:
                element_gaussian_distances = []
                element_distances = []
                for element, psp in self.atomic_psp.items():
                    gaussian_distances = self._get_default_cutoff_element_gaussians(
                        sigma,
                        psp,
                        self.max_num_gaussians,
                        threshold=self.overlap_threshold * factors[order],
                    )
                    element_distance = np.max(gaussian_distances).item()
                    element_distances.append(element_distance)
                    element_gaussian_distances += gaussian_distances
                elemental_order_sigma_gaussian_cutoffs.append(
                    element_gaussian_distances
                )
                elemental_order_sigma_cutoffs.append(element_distances)
                order_sigma_index_dict[(order, sigma)] = order_sigma_index
                order_sigma_index += 1

        self.order_sigma_index_dict = order_sigma_index_dict
        self.nsigmas = len(elemental_order_sigma_gaussian_cutoffs)

        elemental_order_sigma_gaussian_cutoffs = np.asarray(
            elemental_order_sigma_gaussian_cutoffs, dtype=np.float64, order="C"
        )

        elemental_order_sigma_cutoffs = np.asarray(
            elemental_order_sigma_cutoffs, dtype=np.float64, order="C"
        )

        self.params_set[
            "elemental_order_sigma_gaussian_cutoffs"
        ] = elemental_order_sigma_gaussian_cutoffs
        # [[ordersigma1_element1_gaussian1_cutoff, ordersigma1_element1_gaussian2_cutoff,
        #            ...ordersigma1_element2_gaussian1_cutoff,,...]
        #  [ordersigma2_element1_gaussian1_cutoff, ordersigma2_element1_gaussian2_cutoff,
        #            ...ordersigma2_element2_gaussian1_cutoff,,...]
        #  ...
        #  [ordersigman_element1_gaussian1_cutoff, ordersigman_element1_gaussian2_cutoff,
        #            ...ordersigman_element2_gaussian1_cutoff,,...]]
        self.params_set[
            "elemental_order_sigma_gaussian_cutoffs_p"
        ] = _gen_2Darray_for_ffi(elemental_order_sigma_gaussian_cutoffs, ffi)

        self.params_set["elemental_order_sigma_cutoffs"] = elemental_order_sigma_cutoffs
        # [[sigma1_element1_cutoff, sigma1_element2_cutoff,...]
        #  [sigma2_element1_cutoff, sigma2_element2_cutoff,...]
        #  ...
        #  [sigman_element1_cutoff, sigman_element2_cutoff,...]]
        self.params_set["elemental_order_sigma_cutoffs_p"] = _gen_2Darray_for_ffi(
            elemental_order_sigma_cutoffs, ffi
        )
        return

    def _get_sigma_elemental_cutoffs(self):
        """
        private method for getting cutoff distance based on probe gaussian sigma
        and widest elemental gaussian
        """
        elemental_sigma_cutoffs = []
        sigma_index_dict = {}

        for sigma_index, (order, sigma) in enumerate(self.desc_list):
            if sigma in sigma_index_dict:
                continue
            sigma_index_dict[sigma] = sigma_index
            element_distances = []
            for element, psp in self.atomic_psp.items():
                element_distance = self._get_default_cutoff_single_element(
                    sigma, psp, threshold=self.overlap_threshold
                )
                element_distances.append(element_distance)
            elemental_sigma_cutoffs.append(element_distances)

        self.sigma_index_dict = sigma_index_dict
        self.nsigmas = len(elemental_sigma_cutoffs)

        elemental_sigma_cutoffs = np.asarray(
            elemental_sigma_cutoffs, dtype=np.float64, order="C"
        )

        self.params_set["elemental_sigma_cutoffs"] = elemental_sigma_cutoffs
        # [[sigma1_element1_cutoff, sigma1_element2_cutoff,...]
        #  [sigma2_element1_cutoff, sigma2_element2_cutoff,...]
        #  ...
        #  [sigman_element1_cutoff, sigman_element2_cutoff,...]]
        self.params_set["elemental_sigma_cutoffs_p"] = _gen_2Darray_for_ffi(
            elemental_sigma_cutoffs, ffi
        )
        return

    def _get_sigma_cutoffs(self):
        """
        private method for getting cutoff distance based on probe gaussian sigma
        """
        elemental_sigma_cutoffs = []
        result = {}
        for sigma_index, (order, sigma) in enumerate(self.desc_list):
            if sigma in result:
                continue
            element_distances = []
            for element, psp in self.atomic_psp.items():
                element_distance = self._get_default_cutoff_single_element(
                    sigma, psp, threshold=self.overlap_threshold
                )
                element_distances.append(element_distance)
            result[sigma] = np.max(element_distances)

        self.sigma_cutoffs = result
        return

    # *******************************************************
    # setup parameters
    # *******************************************************

    def _get_desc_list(self):
        """
        private method for getting feature setup
        """
        if (
            "GMPs_detailed_list" in self.GMPs
            and self.GMPs["GMPs_detailed_list"] is not None
        ):
            self.desc_list = sorted(self.GMPs["GMPs_detailed_list"])
        else:
            self.desc_list = []
            for order in sorted(self.GMPs["GMPs"]["orders"]):
                for sigma in sorted(self.GMPs["GMPs"]["sigmas"]):
                    if order == -1:
                        self.desc_list.append((-1, 0))
                        break
                    elif order >= 0:
                        self.desc_list.append((order, sigma))
                    else:
                        raise ValueError

    def _load_psp_files(self):
        """
        private method for loading pseudo potential files
        """
        atomic_gaussian_setup = {}
        atomic_psp = {}
        elements = []

        new = False
        temp_params = []
        current_element = None
        with open(self.GMPs["psp_path"], "r") as f:
            Lines = f.readlines()
            for line in Lines:
                line = line.strip()
                if not line.startswith("#"):
                    if line.startswith("!") and line.endswith("elements"):
                        num_elements = int(line.split()[1])
                        continue
                    if line.startswith("*"):
                        new = True
                        if current_element is not None:
                            atomic_gaussian_setup[current_element_index] = np.asarray(
                                temp_params, dtype=np.float64, order="C"
                            )
                            temp_params = []
                        continue
                    if new:
                        current_element = line.split()[0]
                        current_element_index = ATOM_SYMBOL_TO_INDEX_DICT[
                            current_element
                        ]
                        elements.append(current_element)
                        atomic_psp[current_element] = []
                        new = False
                        continue
                    temp = line.split()
                    temp_params += [float(temp[0]), float(temp[1])]
                    atomic_psp[current_element].append([float(temp[0]), float(temp[1])])
            atomic_gaussian_setup[current_element_index] = np.asarray(
                temp_params, dtype=np.float64, order="C"
            )

        self.atomic_gaussian_setup = atomic_gaussian_setup
        self.atomic_psp = atomic_psp

        self.elements = elements
        self.element_indices = list_symbols_to_indices(self.elements)

        max_gaussian_count = 0
        ngaussian_list = list()
        self.params_set = dict()
        for element_index in self.element_indices:
            self.params_set[element_index] = dict()
            self.params_set[element_index][
                "gaussian_params"
            ] = self.atomic_gaussian_setup[element_index]
            self.params_set[element_index]["gaussian_count"] = int(
                len(self.atomic_gaussian_setup[element_index]) / 2
            )
            ngaussian_list.append(self.params_set[element_index]["gaussian_count"])
            # print("self.params_set[element_index]: {}".format(self.params_set[element_index]))

        ngaussian_list = np.asarray(ngaussian_list, dtype=np.intc, order="C")
        max_gaussian_count = np.max(ngaussian_list)
        overall_gaussian_params = list()
        for element_index in self.element_indices:
            temp = np.zeros(max_gaussian_count * 2)
            temp[
                : self.params_set[element_index]["gaussian_count"] * 2
            ] = self.params_set[element_index]["gaussian_params"]
            overall_gaussian_params.append(temp)

        element_index_to_order_list = np.zeros(120, dtype=np.intc)
        for i, element_index in enumerate(self.element_indices):
            element_index_to_order_list[element_index] = i

        overall_gaussian_params = np.asarray(
            overall_gaussian_params, dtype=np.float64, order="C"
        )
        self.params_set["ngaussians"] = ngaussian_list
        self.params_set["ngaussians_p"] = ffi.cast("int *", ngaussian_list.ctypes.data)
        self.params_set["gaussian_params"] = overall_gaussian_params
        self.params_set["gaussian_params_p"] = _gen_2Darray_for_ffi(
            overall_gaussian_params, ffi
        )
        self.params_set["element_index_to_order"] = element_index_to_order_list
        self.params_set["element_index_to_order_p"] = ffi.cast(
            "int *", element_index_to_order_list.ctypes.data
        )

        return

    def _get_feature_setup(self):
        """
        private method for preparing the parameters of feature setup
        """
        self.solid_harmonic = self.GMPs.get("solid_harmonics", True)
        solid_harmonic_i = 1 if self.solid_harmonic else 0
        square = self.GMPs.get("square", False)
        square_i = 1 if square else 0
        rs_setup = self.GMPs.get("rs_setup", {"setup": "constant", "rs": 1.0})
        if rs_setup["setup"] == "scale":
            rs_scale = rs_setup["scale_factor"]
        elif rs_setup["setup"] == "constant":
            rs_scale = -1.0 * rs_setup["rs"]
        else:
            rs_scale = -1.0

        feature_setup = []
        if self.custom_cutoff == -1:
            for order, sigma in self.desc_list:
                if order == -1:
                    feature_setup.append(
                        [
                            order,
                            square_i,
                            solid_harmonic_i,
                            0.0,
                            1.0,
                            1.0,  # A
                            0.0,  # alpha
                            self.cutoff,
                            1.0,  # inv_Rs
                        ]
                    )
                else:
                    feature_setup.append(
                        [
                            order,
                            square_i,
                            solid_harmonic_i,
                            sigma,
                            1.0,
                            (1.0 / (sigma * np.sqrt(2.0 * np.pi))) ** 3,
                            1.0 / (2.0 * sigma * sigma),
                            self.cutoff,
                            (1.0 / (rs_scale * sigma))
                            if rs_scale > 0
                            else (1.0 / (-1.0 * rs_scale)),
                        ]
                    )

        elif self.custom_cutoff == 0:
            for order, sigma in self.desc_list:
                if order == -1:
                    feature_setup.append(
                        [
                            order,
                            square_i,
                            solid_harmonic_i,
                            0.0,
                            1.0,
                            1.0,  # A
                            0.0,  # alpha
                            self.sigma_cutoffs[0],
                            1.0,  # inv_Rs
                        ]
                    )
                else:
                    feature_setup.append(
                        [
                            order,
                            square_i,
                            solid_harmonic_i,
                            sigma,
                            1.0,
                            (1.0 / (sigma * np.sqrt(2.0 * np.pi))) ** 3,
                            1.0 / (2.0 * sigma * sigma),
                            self.sigma_cutoffs[sigma],
                            (1.0 / (rs_scale * sigma))
                            if rs_scale > 0
                            else (1.0 / (-1.0 * rs_scale)),
                        ]
                    )

        elif self.custom_cutoff == 1:
            for order, sigma in self.desc_list:
                if order == -1:
                    feature_setup.append(
                        [
                            order,
                            square_i,
                            solid_harmonic_i,
                            0.0,
                            1.0,
                            1.0,  # A
                            0.0,  # alpha
                            self.sigma_cutoffs[0],
                            1.0,  # inv_Rs
                        ]
                    )
                else:
                    feature_setup.append(
                        [
                            order,
                            square_i,
                            solid_harmonic_i,
                            sigma,
                            1.0,
                            (1.0 / (sigma * np.sqrt(2.0 * np.pi))) ** 3,
                            1.0 / (2.0 * sigma * sigma),
                            self.sigma_cutoffs[sigma],
                            (1.0 / (rs_scale * sigma))
                            if rs_scale > 0
                            else (1.0 / (-1.0 * rs_scale)),
                        ]
                    )

        elif self.custom_cutoff == 2:
            for order, sigma in self.desc_list:
                if order == -1:
                    feature_setup.append(
                        [
                            order,
                            square_i,
                            solid_harmonic_i,
                            self.sigma_index_dict[0],
                            0.0,
                            1.0,
                            1.0,  # A
                            0.0,  # alpha
                            1.0,  # inv_Rs
                        ]
                    )
                else:
                    feature_setup.append(
                        [
                            order,
                            square_i,
                            solid_harmonic_i,
                            self.sigma_index_dict[sigma],
                            sigma,
                            1.0,
                            (1.0 / (sigma * np.sqrt(2.0 * np.pi))) ** 3,
                            1.0 / (2.0 * sigma * sigma),
                            (1.0 / (rs_scale * sigma))
                            if rs_scale > 0
                            else (1.0 / (-1.0 * rs_scale)),
                        ]
                    )

        elif self.custom_cutoff == 3:
            for order, sigma in self.desc_list:
                if order == -1:
                    feature_setup.append(
                        [
                            order,
                            square_i,
                            solid_harmonic_i,
                            self.order_sigma_index_dict[0],
                            0.0,
                            1.0,
                            1.0,  # A
                            0.0,  # alpha
                            1.0,  # inv_Rs
                        ]
                    )
                else:
                    feature_setup.append(
                        [
                            order,
                            square_i,
                            solid_harmonic_i,
                            self.order_sigma_index_dict[(order, sigma)],
                            sigma,
                            1.0,
                            (1.0 / (sigma * np.sqrt(2.0 * np.pi))) ** 3,
                            1.0 / (2.0 * sigma * sigma),
                            (1.0 / (rs_scale * sigma))
                            if rs_scale > 0
                            else (1.0 / (-1.0 * rs_scale)),
                        ]
                    )

        elif self.custom_cutoff == 4:
            for order, sigma in self.desc_list:
                if order == -1:
                    feature_setup.append(
                        [
                            order,
                            square_i,
                            solid_harmonic_i,
                            self.order_sigma_index_dict[0],
                            0.0,
                            1.0,
                            1.0,  # A
                            0.0,  # alpha
                            1.0,  # inv_Rs
                        ]
                    )
                else:
                    feature_setup.append(
                        [
                            order,
                            square_i,
                            solid_harmonic_i,
                            self.order_sigma_index_dict[(order, sigma)],
                            sigma,
                            1.0,
                            (1.0 / (sigma * np.sqrt(2.0 * np.pi))) ** 3,
                            1.0 / (2.0 * sigma * sigma),
                            (1.0 / (rs_scale * sigma))
                            if rs_scale > 0
                            else (1.0 / (-1.0 * rs_scale)),
                        ]
                    )

        else:
            raise NotImplementedError

        self.feature_setup = np.array(feature_setup)

    def _prepare_feature_parameters(self):
        """
        private method for preparing feature parameters
        """
        if self.custom_cutoff <= 1:
            params_i = np.asarray(
                self.feature_setup[:, :3].copy(), dtype=np.intc, order="C"
            )
            params_d = np.asarray(
                self.feature_setup[:, 3:].copy(), dtype=np.float64, order="C"
            )
        elif (
            self.custom_cutoff == 2
            or self.custom_cutoff == 3
            or self.custom_cutoff == 4
        ):
            params_i = np.asarray(
                self.feature_setup[:, :4].copy(), dtype=np.intc, order="C"
            )
            params_d = np.asarray(
                self.feature_setup[:, 4:].copy(), dtype=np.float64, order="C"
            )
        # print(params_d)
        self.params_set["i"] = params_i
        self.params_set["d"] = params_d

        self.params_set["ip"] = _gen_2Darray_for_ffi(self.params_set["i"], ffi, "int")
        self.params_set["dp"] = _gen_2Darray_for_ffi(self.params_set["d"], ffi)
        self.params_set["total"] = np.concatenate(
            (self.params_set["i"], self.params_set["d"]), axis=1
        )
        self.params_set["num"] = len(self.params_set["total"])

        # if "prime_threshold" in self.GMPs:
        #     self.params_set["prime_threshold"] = float(self.GMPs["prime_threshold"])

        self.params_set["log"] = self.GMPs.get("log", False)

        return

    def get_feature_setup_hash(self):
        """
        get hash based on feature setup
        """
        string = ""
        for desc in self.feature_setup:
            for num in desc:
                string += "%.15f" % num
        md5 = hashlib.md5(string.encode("utf-8"))
        hash_result = md5.hexdigest()
        self.feature_setup_hash = hash_result

    def save_feature_setup(self, filename):
        """
        save feature setup to file
        """
        with open(filename, "w") as out_file:
            for desc in self.feature_setup:
                out_file.write("\t".join(str(desc[i]) for i in range(len(desc))) + "\n")

    # *******************************************************
    # calculate features
    # *******************************************************
    def calculate_features(
        self,
        cell,
        pbc,
        atom_positions,
        atom_symbols,
        occupancies,
        ref_positions,
        calc_derivatives,
        calc_occ_derivatives,
    ):
        """
        Method for calculating GMP features based on given image object information
        """
        assert isinstance(cell, np.ndarray)
        assert cell.shape == (3, 3)
        assert isinstance(occupancies, np.ndarray)
        assert isinstance(atom_positions, np.ndarray)
        assert isinstance(ref_positions, np.ndarray)
        assert len(atom_positions) == len(atom_symbols)
        assert len(atom_symbols) == len(occupancies)
        symbols = atom_symbols
        atom_num = len(symbols)
        atom_indices = list_symbols_to_indices(symbols)
        for element_idx in atom_indices:
            if element_idx not in self.element_indices:
                raise ValueError(
                    "****ERROR: contain atoms that don't have pseudo potential"
                    "Please use a different set of pseudo potentials"
                )
        # cell = atoms.cell
        # scaled_ref_positions = cell.scaled_positions(ref_positions)
        scaled_ref_positions = get_scaled_position(cell, ref_positions)
        for number in scaled_ref_positions.flatten():
            if number > 1.0 or number < 0.0:
                raise ValueError(
                    "****ERROR: scaled position not strictly between [0, 1]"
                    "Please check atom position and system cell size are set up correctly"
                )
        scaled_ref_positions = np.array(
            [np.array(v, dtype="float64") for v in scaled_ref_positions],
            dtype="float64",
        )

        atom_indices_p = ffi.cast("int *", atom_indices.ctypes.data)

        occupancies_p = ffi.cast("double *", occupancies.ctypes.data)

        # cart = np.copy(atoms.get_positions(wrap=True), order="C")
        cart = np.copy(atom_positions, order="C")
        # scale = np.copy(atoms.get_scaled_positions(wrap=True), order="C")
        scale = get_scaled_position(cell, atom_positions)
        scale = np.array(
            [np.array(v, dtype="float64") for v in scale],
            dtype="float64",
        )
        atom_indices_p = ffi.cast("int *", atom_indices.ctypes.data)
        cell = np.copy(cell, order="C")
        # pbc = np.copy(atoms.get_pbc()).astype(np.intc)
        pbc = np.copy(pbc).astype(np.intc)

        cart_p = _gen_2Darray_for_ffi(cart, ffi)
        scale_p = _gen_2Darray_for_ffi(scale, ffi)
        cell_p = _gen_2Darray_for_ffi(cell, ffi)
        pbc_p = ffi.cast("int *", pbc.ctypes.data)

        ref_cart_p = _gen_2Darray_for_ffi(ref_positions, ffi)
        ref_scale_p = _gen_2Darray_for_ffi(scaled_ref_positions, ffi)
        cal_num = len(ref_positions)

        size_info = np.array([atom_num, cal_num, self.params_set["num"]])

        if self.solid_harmonic:
            if calc_derivatives is False and calc_occ_derivatives is False:
                x = np.zeros(
                    [cal_num, self.params_set["num"]], dtype=np.float64, order="C"
                )
                x_p = _gen_2Darray_for_ffi(x, ffi)

                if self.custom_cutoff == -1:
                    errno = lib.calculate_solid_gmpordernorm_noderiv_ref(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                    )

                elif self.custom_cutoff == 0:
                    errno = lib.calculate_solid_gmpordernorm_noderiv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                    )

                elif self.custom_cutoff == 1:
                    errno = lib.calculate_solid_gmpordernorm_sigma_cutoff_noderiv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                    )

                elif self.custom_cutoff == 2:
                    errno = (
                        lib.calculate_solid_gmpordernorm_elemental_sigma_cutoff_noderiv(
                            cell_p,
                            cart_p,
                            occupancies_p,
                            ref_cart_p,
                            scale_p,
                            ref_scale_p,
                            pbc_p,
                            atom_indices_p,
                            atom_num,
                            cal_num,
                            self.nsigmas,
                            self.params_set["ip"],
                            self.params_set["dp"],
                            self.params_set["num"],
                            self.params_set["gaussian_params_p"],
                            self.params_set["ngaussians_p"],
                            self.params_set["elemental_sigma_cutoffs_p"],
                            self.params_set["element_index_to_order_p"],
                            x_p,
                        )
                    )

                elif self.custom_cutoff == 3:
                    errno = (
                        lib.calculate_solid_gmpordernorm_elemental_sigma_cutoff_noderiv(
                            cell_p,
                            cart_p,
                            occupancies_p,
                            ref_cart_p,
                            scale_p,
                            ref_scale_p,
                            pbc_p,
                            atom_indices_p,
                            atom_num,
                            cal_num,
                            self.nsigmas,
                            self.params_set["ip"],
                            self.params_set["dp"],
                            self.params_set["num"],
                            self.params_set["gaussian_params_p"],
                            self.params_set["ngaussians_p"],
                            self.params_set["elemental_order_sigma_cutoffs_p"],
                            self.params_set["element_index_to_order_p"],
                            x_p,
                        )
                    )

                elif self.custom_cutoff == 4:
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.nsigmas,
                        self.max_num_gaussians,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["elemental_order_sigma_cutoffs_p"],
                        self.params_set["elemental_order_sigma_gaussian_cutoffs_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                    )

                else:
                    raise NotImplementedError

                if errno == 1:
                    raise NotImplementedError("Feature not implemented!")
                fp = np.array(x, dtype=np.float64)

                return size_info, fp, None, None, None, None, None

            elif calc_derivatives is False and calc_occ_derivatives is True:
                x = np.zeros(
                    [cal_num, self.params_set["num"]], dtype=np.float64, order="C"
                )
                dxdocc = np.zeros(
                    [cal_num * self.params_set["num"], atom_num],
                    dtype=np.float64,
                    order="C",
                )

                x_p = _gen_2Darray_for_ffi(x, ffi)
                dxdocc_p = _gen_2Darray_for_ffi(dxdocc, ffi)

                if self.custom_cutoff == -1:
                    errno = lib.calculate_solid_gmpordernorm_occ_deriv_ref(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dxdocc_p,
                    )

                elif self.custom_cutoff == 0:
                    errno = lib.calculate_solid_gmpordernorm_occ_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dxdocc_p,
                    )

                elif self.custom_cutoff == 1:
                    errno = lib.calculate_solid_gmpordernorm_sigma_cutoff_occ_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dxdocc_p,
                    )

                elif self.custom_cutoff == 2:
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_cutoff_occ_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.nsigmas,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["elemental_sigma_cutoffs_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dxdocc_p,
                    )

                elif self.custom_cutoff == 3:
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_cutoff_occ_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.nsigmas,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["elemental_order_sigma_cutoffs_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dxdocc_p,
                    )

                elif self.custom_cutoff == 4:
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_occ_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.nsigmas,
                        self.max_num_gaussians,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["elemental_order_sigma_cutoffs_p"],
                        self.params_set["elemental_order_sigma_gaussian_cutoffs_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dxdocc_p,
                    )

                else:
                    raise NotImplementedError
                if errno == 1:
                    raise NotImplementedError("Feature not implemented!")
                fp = np.array(x, dtype=np.float64)
                dfp_docc = np.array(dxdocc, dtype=np.float64)

                return size_info, fp, None, None, None, None, dfp_docc

            elif calc_derivatives is True and calc_occ_derivatives is False:
                x = np.zeros(
                    [cal_num, self.params_set["num"]], dtype=np.float64, order="C"
                )
                dx = np.zeros(
                    [cal_num * self.params_set["num"], atom_num * 3],
                    dtype=np.float64,
                    order="C",
                )

                x_p = _gen_2Darray_for_ffi(x, ffi)
                dx_p = _gen_2Darray_for_ffi(dx, ffi)

                if self.custom_cutoff == -1:
                    errno = lib.calculate_solid_gmpordernorm_fp_deriv_ref(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                    )

                elif self.custom_cutoff == 0:
                    errno = lib.calculate_solid_gmpordernorm_fp_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                    )

                elif self.custom_cutoff == 1:
                    errno = lib.calculate_solid_gmpordernorm_sigma_cutoff_fp_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                    )

                elif self.custom_cutoff == 2:
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_cutoff_fp_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.nsigmas,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["elemental_sigma_cutoffs_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                    )

                elif self.custom_cutoff == 3:
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_cutoff_fp_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.nsigmas,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["elemental_order_sigma_cutoffs_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                    )

                elif self.custom_cutoff == 4:
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_fp_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.nsigmas,
                        self.max_num_gaussians,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["elemental_order_sigma_cutoffs_p"],
                        self.params_set["elemental_order_sigma_gaussian_cutoffs_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                    )

                else:
                    raise NotImplementedError

                if errno == 1:
                    raise NotImplementedError("Feature not implemented!")
                fp = np.array(x, dtype=np.float64)
                fp_prime = np.array(dx, dtype=np.float64)

                scipy_sparse_fp_prime = sparse.coo_matrix(fp_prime)

                return (
                    size_info,
                    fp,
                    scipy_sparse_fp_prime.data,
                    scipy_sparse_fp_prime.row,
                    scipy_sparse_fp_prime.col,
                    np.array(fp_prime.shape),
                    None,
                )

            else:
                x = np.zeros(
                    [cal_num, self.params_set["num"]], dtype=np.float64, order="C"
                )
                dx = np.zeros(
                    [cal_num * self.params_set["num"], atom_num * 3],
                    dtype=np.float64,
                    order="C",
                )
                dxdocc = np.zeros(
                    [cal_num * self.params_set["num"], atom_num],
                    dtype=np.float64,
                    order="C",
                )

                x_p = _gen_2Darray_for_ffi(x, ffi)
                dx_p = _gen_2Darray_for_ffi(dx, ffi)
                dxdocc_p = _gen_2Darray_for_ffi(dxdocc, ffi)

                if self.custom_cutoff == -1:
                    errno = lib.calculate_solid_gmpordernorm_fp_occ_deriv_ref(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                        dxdocc_p,
                    )

                elif self.custom_cutoff == 0:
                    errno = lib.calculate_solid_gmpordernorm_fp_occ_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                        dxdocc_p,
                    )

                elif self.custom_cutoff == 1:
                    errno = lib.calculate_solid_gmpordernorm_sigma_cutoff_fp_occ_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                        dxdocc_p,
                    )

                elif self.custom_cutoff == 2:
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_cutoff_fp_occ_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.nsigmas,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["elemental_sigma_cutoffs_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                        dxdocc_p,
                    )

                elif self.custom_cutoff == 3:
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_cutoff_fp_occ_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.nsigmas,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["elemental_order_sigma_cutoffs_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                        dxdocc_p,
                    )

                elif self.custom_cutoff == 4:
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_fp_occ_deriv(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.nsigmas,
                        self.max_num_gaussians,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["elemental_order_sigma_cutoffs_p"],
                        self.params_set["elemental_order_sigma_gaussian_cutoffs_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                        dxdocc_p,
                    )

                else:
                    raise NotImplementedError

                if errno == 1:
                    raise NotImplementedError("Feature not implemented!")
                fp = np.array(x, dtype=np.float64)
                fp_prime = np.array(dx, dtype=np.float64)
                dfp_docc = np.array(dxdocc, dtype=np.float64)

                scipy_sparse_fp_prime = sparse.coo_matrix(fp_prime)

                return (
                    size_info,
                    fp,
                    scipy_sparse_fp_prime.data,
                    scipy_sparse_fp_prime.row,
                    scipy_sparse_fp_prime.col,
                    np.array(fp_prime.shape),
                    dfp_docc,
                )

        else:

            if calc_derivatives is False and calc_occ_derivatives is False:
                x = np.zeros(
                    [cal_num, self.params_set["num"]], dtype=np.float64, order="C"
                )
                x_p = _gen_2Darray_for_ffi(x, ffi)
                if self.custom_cutoff == -1:
                    errno = lib.calculate_surface_gmpordernorm_noderiv_ref(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                    )

                else:
                    raise NotImplementedError

                if errno == 1:
                    raise NotImplementedError("Feature not implemented!")
                fp = np.array(x, dtype=np.float64)

                return (
                    size_info,
                    fp,
                    None,
                    None,
                    None,
                    None,
                    None,
                )

            elif calc_derivatives is False and calc_occ_derivatives is True:
                raise NotImplementedError

            elif calc_derivatives is True and calc_occ_derivatives is False:
                x = np.zeros(
                    [cal_num, self.params_set["num"]], dtype=np.float64, order="C"
                )
                dx = np.zeros(
                    [cal_num * self.params_set["num"], atom_num * 3],
                    dtype=np.float64,
                    order="C",
                )

                x_p = _gen_2Darray_for_ffi(x, ffi)
                dx_p = _gen_2Darray_for_ffi(dx, ffi)

                if self.custom_cutoff == -1:
                    errno = lib.calculate_surface_gmpordernorm_fp_deriv_ref(
                        cell_p,
                        cart_p,
                        occupancies_p,
                        ref_cart_p,
                        scale_p,
                        ref_scale_p,
                        pbc_p,
                        atom_indices_p,
                        atom_num,
                        cal_num,
                        self.params_set["ip"],
                        self.params_set["dp"],
                        self.params_set["num"],
                        self.params_set["gaussian_params_p"],
                        self.params_set["ngaussians_p"],
                        self.params_set["element_index_to_order_p"],
                        x_p,
                        dx_p,
                    )

                else:
                    raise NotImplementedError

                if errno == 1:
                    raise NotImplementedError("Feature not implemented!")
                fp = np.array(x, dtype=np.float64)
                fp_prime = np.array(dx, dtype=np.float64)

                scipy_sparse_fp_prime = sparse.coo_matrix(fp_prime)

                return (
                    size_info,
                    fp,
                    scipy_sparse_fp_prime.data,
                    scipy_sparse_fp_prime.row,
                    scipy_sparse_fp_prime.col,
                    np.array(fp_prime.shape),
                    None,
                )

            else:
                raise NotImplementedError
