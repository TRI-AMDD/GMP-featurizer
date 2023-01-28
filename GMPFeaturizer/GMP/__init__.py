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
    def __init__(
        self,
        GMPs,
        feature_database="cache/features/",
    ):
        super().__init__()
        self.feature_type = "GMP"
        self.feature_database = feature_database
        self.GMPs = GMPs
        self.custom_cutoff = self.GMPs.get("custom_cutoff", 4)

        self._load_psp_files()
        self.overlap_threshold = self.GMPs.get("overlap_threshold", 1e-11)
        self._get_cutoffs()
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
            # print(temp_distance)
            distances_list.append(temp_distance)

        return np.max(distances_list).item()

    def _get_sigma_list(self, attach_zero=True):
        if "GMPs_detailed_list" in self.GMPs:
            raise NotImplementedError
        else:
            sigmas = list(self.GMPs["GMPs"]["sigmas"]).copy()
            if attach_zero and -1 in self.GMPs["GMPs"]["orders"]:
                sigmas.append(0)
        return sigmas

    def _get_order_list(self):
        if "GMPs_detailed_list" in self.GMPs:
            raise NotImplementedError
        orders = list(self.GMPs["GMPs"]["orders"])
        return orders

    def _get_cutoffs(self):
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
        # self.params_set["element_index_to_order"]
        elemental_order_sigma_cutoffs = []
        sigmas = self._get_sigma_list(attach_zero=False)
        orders = self._get_order_list()

        # sigma_index_dict = {}
        order_sigma_index_dict = {}
        result = {}
        order_sigma_index = 0
        for order in orders:
            if order == -1:
                element_distances = []
                for element, psp in self.atomic_psp.items():
                    element_distance = self._get_default_cutoff_single_element(
                        0, psp, threshold=self.overlap_threshold
                    )
                    element_distances.append(element_distance)
                    print(order, 0, element, element_distance)
                elemental_order_sigma_cutoffs.append(element_distances)
                order_sigma_index_dict[0] = order_sigma_index
                order_sigma_index += 1

            else:
                for sigma in sigmas:
                    element_distances = []
                    for element, psp in self.atomic_psp.items():
                        element_distance = self._get_default_cutoff_single_element(
                            sigma,
                            psp,
                            threshold=self.overlap_threshold * factors[order],
                        )
                        element_distances.append(element_distance)
                        print(order, sigma, element, element_distance)
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
        # self.params_set["element_index_to_order"]
        elemental_order_sigma_cutoffs = []
        elemental_order_sigma_gaussian_cutoffs = []
        sigmas = self._get_sigma_list(attach_zero=False)
        orders = self._get_order_list()

        self.max_num_gaussians = 0
        for psp in self.atomic_psp.values():
            if len(psp) > self.max_num_gaussians:
                self.max_num_gaussians = len(psp)

        # sigma_index_dict = {}
        order_sigma_index_dict = {}
        result = {}
        order_sigma_index = 0
        for order in orders:
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
                for sigma in sigmas:
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
                        # print(order, sigma, element, element_distances)
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
        # [[ordersigma1_element1_gaussian1_cutoff, ordersigma1_element1_gaussian2_cutoff,...ordersigma1_element2_gaussian1_cutoff,,...]
        #  [ordersigma2_element1_gaussian1_cutoff, ordersigma2_element1_gaussian2_cutoff,...ordersigma2_element2_gaussian1_cutoff,,...]
        #  ...
        #  [ordersigman_element1_gaussian1_cutoff, ordersigman_element1_gaussian2_cutoff,...ordersigman_element2_gaussian1_cutoff,,...]]
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
        elemental_sigma_cutoffs = []
        sigmas = self._get_sigma_list()
        orders = self._get_order_list()
        self.nsigmas = len(sigmas)
        sigma_index_dict = {}

        for sigma_index, sigma in enumerate(sigmas):
            sigma_index_dict[sigma] = sigma_index
            element_distances = []
            for element, psp in self.atomic_psp.items():
                element_distance = self._get_default_cutoff_single_element(
                    sigma, psp, threshold=self.overlap_threshold
                )
                element_distances.append(element_distance)
            elemental_sigma_cutoffs.append(element_distances)

        self.sigma_index_dict = sigma_index_dict
        # print(sigma_index_dict)

        elemental_sigma_cutoffs = np.asarray(
            elemental_sigma_cutoffs, dtype=np.float64, order="C"
        )

        # print(elemental_sigma_cutoffs)

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
        elemental_sigma_cutoffs = []
        sigmas = self._get_sigma_list()
        result = {}
        for sigma_index, sigma in enumerate(sigmas):
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

    def _load_psp_files(self):
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
                    if new == True:
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

        # for element in self.elements:
        #     params = list()
        #     atomic_psp[element] = []
        #     # count = 0
        #     filename = self.GMPs["atom_gaussians"][element]
        #     with open(filename, "r") as fil:
        #         for line in fil:
        #             tmp = line.split()
        #             params += [float(tmp[0]), float(tmp[1])]
        #             atomic_psp[element].append([float(tmp[0]), float(tmp[1])])
        #             # count += 1
        #     element_index = ATOM_SYMBOL_TO_INDEX_DICT[element]
        #     params = np.asarray(params, dtype=np.float64, order="C")
        #     atomic_gaussian_setup[element_index] = params

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

        # print("ngaussian_list: {}".format(ngaussian_list))
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
            # element_index_to_order_list.append(element_index)
            element_index_to_order_list[element_index] = i

        # element_index_to_order_list = np.asarray(element_index_to_order_list, dtype=np.intc, order='C')
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

        # print("ngaussians: {}".format(self.params_set["ngaussians"]))
        # print("gaussian_params: {}".format(self.params_set["gaussian_params"]))
        return

    def _prepare_feature_parameters(self):
        feature_setup = []
        # cutoff = self.GMPs["cutoff"]
        self.solid_harmonic = self.GMPs.get("solid_harmonics", False)
        solid_harmonic_i = 1 if self.solid_harmonic else 0
        square = self.GMPs.get("square", True)
        square_i = 1 if square else 0
        rs_setup = self.GMPs.get("rs_setup", {"setup": "constant", "rs": 1.0})
        if rs_setup["setup"] == "scale":
            rs_scale = rs_setup["scale_factor"]
        elif rs_setup["setup"] == "constant":
            rs_scale = -1.0 * rs_setup["rs"]
        else:
            rs_scale = -1.0

        if "GMPs_detailed_list" in self.GMPs:
            raise NotImplementedError
            # for key in sorted(self.GMPs["GMPs_detailed_list"].keys()):
            #     detail_setup = self.GMPs["GMPs_detailed_list"][key]
            #     sigmas = detail_setup["sigmas"]
            #     feature_setup += [
            #         [
            #             detail_setup["order"],
            #             square_i,
            #             solid_harmonic_i,
            #             sigma,
            #             1.0,
            #             (1.0 / (sigma * np.sqrt(2.0 * np.pi))) ** 3,
            #             1.0 / (2.0 * sigma * sigma),
            #             self.get_cutoff(sigma),
            #             (1.0 / (rs_scale * sigma))
            #             if rs_scale > 0
            #             else (1.0 / (-1.0 * rs_scale)),
            #         ]
            #         for sigma in sigmas
            #     ]

        else:
            if self.custom_cutoff == -1:
                for i in self.GMPs["GMPs"]["orders"]:
                    if i == -1:
                        feature_setup += [
                            [
                                i,
                                square_i,
                                solid_harmonic_i,
                                0.0,
                                1.0,
                                1.0,  # A
                                0.0,  # alpha
                                self.cutoff,
                                1.0,  # inv_Rs
                            ]
                        ]
                    else:
                        feature_setup += [
                            [
                                i,
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
                            for sigma in self.GMPs["GMPs"]["sigmas"]
                        ]

            elif self.custom_cutoff == 0:
                for i in self.GMPs["GMPs"]["orders"]:
                    if i == -1:
                        feature_setup += [
                            [
                                i,
                                square_i,
                                solid_harmonic_i,
                                0.0,
                                1.0,
                                1.0,  # A
                                0.0,  # alpha
                                self.sigma_cutoffs[0],
                                1.0,  # inv_Rs
                            ]
                        ]
                    else:
                        feature_setup += [
                            [
                                i,
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
                            for sigma in self.GMPs["GMPs"]["sigmas"]
                        ]

            elif self.custom_cutoff == 1:
                for i in self.GMPs["GMPs"]["orders"]:
                    if i == -1:
                        feature_setup += [
                            [
                                i,
                                square_i,
                                solid_harmonic_i,
                                0.0,
                                1.0,
                                1.0,  # A
                                0.0,  # alpha
                                self.sigma_cutoffs[0],
                                1.0,  # inv_Rs
                            ]
                        ]
                    else:
                        feature_setup += [
                            [
                                i,
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
                            for sigma in self.GMPs["GMPs"]["sigmas"]
                        ]

            elif self.custom_cutoff == 2:
                for i in self.GMPs["GMPs"]["orders"]:
                    if i == -1:
                        feature_setup += [
                            [
                                i,
                                square_i,
                                solid_harmonic_i,
                                self.sigma_index_dict[0],
                                0.0,
                                1.0,
                                1.0,  # A
                                0.0,  # alpha
                                1.0,  # inv_Rs
                            ]
                        ]
                    else:
                        feature_setup += [
                            [
                                i,
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
                            for sigma in self.GMPs["GMPs"]["sigmas"]
                        ]

            elif self.custom_cutoff == 3:
                for i in self.GMPs["GMPs"]["orders"]:
                    if i == -1:
                        feature_setup += [
                            [
                                i,
                                square_i,
                                solid_harmonic_i,
                                self.order_sigma_index_dict[0],
                                0.0,
                                1.0,
                                1.0,  # A
                                0.0,  # alpha
                                1.0,  # inv_Rs
                            ]
                        ]
                    else:
                        feature_setup += [
                            [
                                i,
                                square_i,
                                solid_harmonic_i,
                                self.order_sigma_index_dict[(i, sigma)],
                                sigma,
                                1.0,
                                (1.0 / (sigma * np.sqrt(2.0 * np.pi))) ** 3,
                                1.0 / (2.0 * sigma * sigma),
                                (1.0 / (rs_scale * sigma))
                                if rs_scale > 0
                                else (1.0 / (-1.0 * rs_scale)),
                            ]
                            for sigma in self.GMPs["GMPs"]["sigmas"]
                        ]
            elif self.custom_cutoff == 4:
                for i in self.GMPs["GMPs"]["orders"]:
                    if i == -1:
                        feature_setup += [
                            [
                                i,
                                square_i,
                                solid_harmonic_i,
                                self.order_sigma_index_dict[0],
                                0.0,
                                1.0,
                                1.0,  # A
                                0.0,  # alpha
                                1.0,  # inv_Rs
                            ]
                        ]
                    else:
                        feature_setup += [
                            [
                                i,
                                square_i,
                                solid_harmonic_i,
                                self.order_sigma_index_dict[(i, sigma)],
                                sigma,
                                1.0,
                                (1.0 / (sigma * np.sqrt(2.0 * np.pi))) ** 3,
                                1.0 / (2.0 * sigma * sigma),
                                (1.0 / (rs_scale * sigma))
                                if rs_scale > 0
                                else (1.0 / (-1.0 * rs_scale)),
                            ]
                            for sigma in self.GMPs["GMPs"]["sigmas"]
                        ]
            else:
                raise NotImplementedError

        self.feature_setup = np.array(feature_setup)
        # print(self.feature_setup)
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

        if "prime_threshold" in self.GMPs:
            self.params_set["prime_threshold"] = float(self.GMPs["prime_threshold"])

        self.params_set["log"] = self.GMPs.get("log", False)

        return

    def get_feature_setup_hash(self):
        # set self.feature_setup_hash
        string = ""
        for desc in self.feature_setup:
            for num in desc:
                string += "%.15f" % num
        md5 = hashlib.md5(string.encode("utf-8"))
        hash_result = md5.hexdigest()
        self.feature_setup_hash = hash_result

    def save_feature_setup(self, filename):
        with open(filename, "w") as out_file:
            for desc in self.feature_setup:
                out_file.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        int(desc[0]),
                        int(desc[1]),
                        desc[2],
                        desc[3],
                        desc[4],
                        desc[5],
                        desc[6],
                    )
                )

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
        calc_derivatives_occ=False,
    ):
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

        # def calculate_features(self, atoms, ref_positions, calc_derivatives):

        #     symbols = np.array(atoms.get_chemical_symbols())
        #     atom_num = len(symbols)
        #     atom_indices = list_symbols_to_indices(symbols)
        #     cell = atoms.cell
        #     scaled_ref_positions = cell.scaled_positions(ref_positions)
        #     scaled_ref_positions = np.array(
        #         [np.array(v, dtype="float64") for v in scaled_ref_positions],
        #         dtype="float64",
        #     )
        #     atom_indices_p = ffi.cast("int *", atom_indices.ctypes.data)

        #     cart = np.copy(atoms.get_positions(wrap=True), order="C")
        #     scale = np.copy(atoms.get_scaled_positions(wrap=True), order="C")
        #     cell = np.copy(atoms.cell, order="C")
        #     pbc = np.copy(atoms.get_pbc()).astype(np.intc)

        #     cart_p = _gen_2Darray_for_ffi(cart, ffi)
        #     scale_p = _gen_2Darray_for_ffi(scale, ffi)
        #     cell_p = _gen_2Darray_for_ffi(cell, ffi)
        #     pbc_p = ffi.cast("int *", pbc.ctypes.data)

        #     ref_cart_p = _gen_2Darray_for_ffi(ref_positions, ffi)
        #     ref_scale_p = _gen_2Darray_for_ffi(scaled_ref_positions, ffi)
        #     cal_num = len(ref_positions)

        #     size_info = np.array([atom_num, cal_num, self.params_set["num"]])

        if calc_derivatives_occ:
            x = np.zeros([cal_num, self.params_set["num"]], dtype=np.float64, order="C")
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

            if self.solid_harmonic:
                if self.custom_cutoff == 3:
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_cutoff_occ_derivative(
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
                else:
                    raise NotImplementedError

            else:
                raise NotImplementedError

            if errno == 1:
                raise NotImplementedError("Feature not implemented!")
            fp = np.array(x, dtype=np.float64)
            fp_prime = np.array(dx, dtype=np.float64)
            fp_occ_prime = np.array(dxdocc, dtype=np.float64)

            scipy_sparse_fp_prime = sparse.coo_matrix(fp_prime)

            return (
                size_info,
                fp,
                scipy_sparse_fp_prime.data,
                scipy_sparse_fp_prime.row,
                scipy_sparse_fp_prime.col,
                np.array(fp_prime.shape),
                fp_occ_prime,
            )

        if calc_derivatives:

            x = np.zeros([cal_num, self.params_set["num"]], dtype=np.float64, order="C")
            dx = np.zeros(
                [cal_num * self.params_set["num"], atom_num * 3],
                dtype=np.float64,
                order="C",
            )

            x_p = _gen_2Darray_for_ffi(x, ffi)
            dx_p = _gen_2Darray_for_ffi(dx, ffi)

            if self.solid_harmonic:
                if self.custom_cutoff == 0:
                    errno = lib.calculate_solid_gmpordernorm(
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
                elif self.custom_cutoff == 3:
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_cutoff(
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
                    errno = lib.calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff(
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
                errno = lib.calculate_gmpordernorm(
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

            if errno == 1:
                raise NotImplementedError("Feature not implemented!")
            fp = np.array(x, dtype=np.float64)
            fp_prime = np.array(dx, dtype=np.float64)

            # threshold = 1e-9
            # super_threshold_indices_prime = np.abs(fp_prime) < threshold
            # print("threshold: {} \tnum points set to zero:{} \t outof: {}".format(threshold, np.sum(super_threshold_indices_prime), fp_prime.shape[0] * fp_prime.shape[1]))
            # fp_prime[super_threshold_indices_prime] = 0.0
            # print(fp_prime)
            scipy_sparse_fp_prime = sparse.coo_matrix(fp_prime)
            # print(fp)
            # print(fp.shape)
            # print(np.sum(super_threshold_indices))
            # print(np.min(np.abs(scipy_sparse_fp_prime.data)))
            # print("density: {}% \n\n----------------------".format(100*len(scipy_sparse_fp_prime.data) / (fp_prime.shape[0] * fp_prime.shape[1])))
            # if self.params_set["log"]:
            #     raise NotImplementedError

            return (
                size_info,
                fp,
                scipy_sparse_fp_prime.data,
                scipy_sparse_fp_prime.row,
                scipy_sparse_fp_prime.col,
                np.array(fp_prime.shape),
            )

        else:
            x = np.zeros([cal_num, self.params_set["num"]], dtype=np.float64, order="C")

            x_p = _gen_2Darray_for_ffi(x, ffi)

            if self.solid_harmonic:
                if self.custom_cutoff == 1:
                    errno = lib.calculate_solid_gmpordernorm_noderiv_sigma_cutoff(
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
                        lib.calculate_solid_gmpordernorm_noderiv_elemental_sigma_cutoff(
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
                        lib.calculate_solid_gmpordernorm_noderiv_elemental_sigma_cutoff(
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
                    errno = lib.calculate_solid_gmpordernorm_noderiv_elemental_sigma_gaussian_cutoff(
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
                elif self.custom_cutoff == 0:
                    errno = lib.calculate_solid_gmpordernorm_noderiv_original(
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
                elif self.custom_cutoff == -1:
                    errno = lib.calculate_solid_gmpordernorm_noderiv_original(
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

            else:
                errno = lib.calculate_gmpordernorm_noderiv(
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

            if errno == 1:
                raise NotImplementedError("Feature not implemented!")

            # print("calculation done")

            fp = np.array(x, dtype=np.float64)
            # if self.params_set["log"]:
            #     fp = np.abs(fp)
            #     fp[fp < 1e-8] = 1e-8
            #     fp = np.log10(fp)

            return size_info, fp, None, None, None, None
