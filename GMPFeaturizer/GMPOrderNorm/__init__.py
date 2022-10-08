import hashlib
import math

import numpy as np
from scipy import sparse

from ..base_feature import BaseFeature
from ..constants import ATOM_SYMBOL_TO_INDEX_DICT
from ..util import _gen_2Darray_for_ffi, list_symbols_to_indices, _gen_2Darray_for_ffi2
from ._libgmpordernorm import ffi, lib


class GMPOrderNorm(BaseFeature):
    def __init__(
        self,
        GMPs,
        elements,
        # mode="atom-centered",
    ):
        super().__init__()
        self.feature_type = "GMPOrderNorm"
        self.GMPs = GMPs
        self.custom_cutoff = self.GMPs.get("custom_cutoff",2)
        self.elements = elements
        self.element_indices = list_symbols_to_indices(elements)

        self._load_psp_files()
        self.overlap_threshold = self.GMPs.get("overlap_threshold", 1e-6)
        self._get_sigma_cutoffs()
        self._prepare_feature_parameters()

        
        # if "cutoff" not in self.GMPs:
        #     self.set_max_cutoff()

        # self.prepare_feature_parameters()

        # self.get_feature_setup_hash()

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

    def _get_C1_C2(self,A,B,alpha,beta):
        if alpha == 0:
            C1 = B
            C2 = -1.0 * beta
        else:
            temp = math.sqrt(math.pi / (alpha + beta))
            C1 = A * B * temp * temp * temp
            C2 = -1.0 * (alpha * beta / (alpha + beta))
        return C1, C2
    
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
            B,beta = gaussian
            C1, C2 = self._get_C1_C2(A,B,alpha,beta)
            C1 = abs(C1)
            temp_distance = math.sqrt((log_threshold - C1) / C2)
            # print(temp_distance)
            distances_list.append(temp_distance)

        return np.max(distances_list).item()
    
    def _get_sigma_list(self):
        if "GMPs_detailed_list" in self.GMPs:
            raise NotImplementedError
        else:
            sigmas = list(self.GMPs["GMPs"]["sigmas"]).copy()
            if -1 in self.GMPs["GMPs"]["orders"]:
                sigmas.append(0)
        return sigmas 

    def _get_sigma_cutoffs(self):
        self.params_set["element_index_to_order"]
        elemental_sigma_cutoffs = []
        sigmas = self._get_sigma_list()
        self.nsigmas = len(sigmas)
        sigma_index_dict = {}
        result = {}
        for sigma_index, sigma in enumerate(sigmas):
            sigma_index_dict[sigma] = sigma_index
            element_distances = []
            for element,psp in self.atomic_psp.items():
                element_distance = self._get_default_cutoff_single_element(sigma, psp, threshold=self.overlap_threshold)
                element_distances.append(element_distance)
                # print(sigma, element, element_distance)
            elemental_sigma_cutoffs.append(element_distances)
            result[sigma] = np.max(element_distances)
            
        self.sigma_cutoffs = result
        # print(self.sigma_cutoffs)
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

    # def set_max_cutoff(self, threshold=1e-8):
    #     max_sigma = 0.0
    #     if "GMPs_detailed_list" in self.GMPs:
    #         for detail_setup in self.GMPs["GMPs_detailed_list"]:
    #             sigmas = detail_setup["sigmas"]
    #             max_sigma = max(max(sigmas), max_sigma)

    #     else:
    #         sigmas = self.GMPs["GMPs"]["sigmas"]
    #         max_sigma = max(max(sigmas), max_sigma)

    #     max_cutoff = self.get_cutoff(max_sigma)

    #     print("max cutoff distance: {} A".format(max_cutoff))

    #     # self.GMPs["max_cutoff"] = max_cutoff
    #     self.max_cutoff = max_cutoff

    #     return
    
    # def get_cutoff(self, sigma):

    #     A = 1.0 / (sigma * np.sqrt(2.0 * math.pi))
    #     alpha = 1.0 / (2.0 * sigma * sigma)

    #     cutoff = round(np.sqrt(np.log(self.cutoff_threshold / A) / (-alpha)), 3)

    # def _get_single_overlap(A,B,alpha,beta):
    #     temp = math.sqrt(M_PI / (alpha + beta))
    #     C1 = A * B * temp * temp * temp
    #     C2 = -1.0 * (alpha * beta / (alpha + beta))
    #     return C1 * exp( C2 * r0_sqr)


    # *******************************************************
    # setup parameters
    # *******************************************************

    def _load_psp_files(self):
        atomic_gaussian_setup = {}
        atomic_psp = {}
        for element in self.elements:
            params = list()
            atomic_psp[element] = []
            # count = 0
            filename = self.GMPs["atom_gaussians"][element]
            with open(filename, "r") as fil:
                for line in fil:
                    tmp = line.split()
                    params += [float(tmp[0]), float(tmp[1])]
                    atomic_psp[element].append([float(tmp[0]), float(tmp[1])])
                    # count += 1
            element_index = ATOM_SYMBOL_TO_INDEX_DICT[element]
            params = np.asarray(params, dtype=np.float64, order="C")
            atomic_gaussian_setup[element_index] = params

        self.atomic_gaussian_setup = atomic_gaussian_setup
        self.atomic_psp = atomic_psp

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
            if self.custom_cutoff <= 1:
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
        elif self.custom_cutoff == 2:
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

    def calculate_features(self, atoms, ref_positions, calc_derivatives):

        symbols = np.array(atoms.get_chemical_symbols())
        atom_num = len(symbols)
        atom_indices = list_symbols_to_indices(symbols)
        cell = atoms.cell
        scaled_ref_positions = cell.scaled_positions(ref_positions)
        scaled_ref_positions = np.array(
            [np.array(v, dtype="float64") for v in scaled_ref_positions],
            dtype="float64",
        )
        atom_indices_p = ffi.cast("int *", atom_indices.ctypes.data)

        cart = np.copy(atoms.get_positions(wrap=True), order="C")
        scale = np.copy(atoms.get_scaled_positions(wrap=True), order="C")
        cell = np.copy(atoms.cell, order="C")
        pbc = np.copy(atoms.get_pbc()).astype(np.intc)

        cart_p = _gen_2Darray_for_ffi(cart, ffi)
        scale_p = _gen_2Darray_for_ffi(scale, ffi)
        cell_p = _gen_2Darray_for_ffi(cell, ffi)
        pbc_p = ffi.cast("int *", pbc.ctypes.data)

        ref_cart_p = _gen_2Darray_for_ffi(ref_positions, ffi)
        ref_scale_p = _gen_2Darray_for_ffi(scaled_ref_positions, ffi)
        cal_num = len(ref_positions)

        size_info = np.array([atom_num, cal_num, self.params_set["num"]])

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
                errno = lib.calculate_solid_gmpordernorm(
                    cell_p,
                    cart_p,
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
                errno = lib.calculate_gmpordernorm(
                    cell_p,
                    cart_p,
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
                    errno = lib.calculate_solid_gmpordernorm_noderiv_elemental_sigma_cutoff(
                        cell_p,
                        cart_p,
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
                elif self.custom_cutoff == 0:
                    errno = lib.calculate_solid_gmpordernorm_noderiv_original(
                        cell_p,
                        cart_p,
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
