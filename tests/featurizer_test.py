import unittest
import numpy as np
from GMPFeaturizer import GMPFeaturizer
import pickle


class FeatruizerTest(unittest.TestCase):
    def base_molecule_featurizer_test(self, GMPs, cores):
        # ensuring the reference method give expected outputs
        with open("./test_files/molecules.p", "rb") as f:
            images = pickle.load(f)

        with open("./test_files/molecules_ref_features.p", "rb") as f:
            reference_features = pickle.load(f)

        featurizer_noderiv = GMPFeaturizer(
            GMPs=GMPs,
            calc_derivatives=False,
            calc_occ_derivatives=False,
            verbose=False,
        )
        features_noderiv = featurizer_noderiv.prepare_features(images, cores=cores)

        featurizer_fp_deriv = GMPFeaturizer(
            GMPs=GMPs,
            calc_derivatives=True,
            calc_occ_derivatives=False,
            verbose=False,
        )
        features_fp_deriv = featurizer_fp_deriv.prepare_features(images, cores=cores)

        featurizer_occ_deriv = GMPFeaturizer(
            GMPs=GMPs,
            calc_derivatives=False,
            calc_occ_derivatives=True,
            verbose=False,
        )
        features_occ_deriv = featurizer_occ_deriv.prepare_features(images, cores=cores)

        featurizer_fp_occ_deriv = GMPFeaturizer(
            GMPs=GMPs,
            calc_derivatives=True,
            calc_occ_derivatives=True,
            verbose=False,
        )
        features_fp_occ_deriv = featurizer_fp_occ_deriv.prepare_features(
            images, cores=cores
        )

        assert np.all(
            [
                np.allclose(
                    features_noderiv[i]["features"],
                    reference_features[i]["features"],
                    rtol=1e-10,
                    atol=1e-15,
                )
                for i in range(len(images))
            ]
        )
        assert np.all(
            [
                np.allclose(
                    features_fp_deriv[i]["features"],
                    reference_features[i]["features"],
                    rtol=1e-10,
                    atol=1e-15,
                )
                for i in range(len(images))
            ]
        )
        assert np.all(
            [
                np.allclose(
                    features_occ_deriv[i]["features"],
                    reference_features[i]["features"],
                    rtol=1e-10,
                    atol=1e-15,
                )
                for i in range(len(images))
            ]
        )
        assert np.all(
            [
                np.allclose(
                    features_fp_occ_deriv[i]["features"],
                    reference_features[i]["features"],
                    rtol=1e-10,
                    atol=1e-15,
                )
                for i in range(len(images))
            ]
        )

    def test_features_ref(self):

        GMPs = {
            "GMPs": {
                "orders": [-1, 0, 1, 2, 3],
                "sigmas": [0.2, 0.4, 0.6, 0.8, 1.0],
            },
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": -1,
            "cutoff": 60,
        }

        self.base_molecule_featurizer_test(GMPs, 1)

    def test_features_method_0(self):

        GMPs = {
            "GMPs": {
                "orders": [-1, 0, 1, 2, 3],
                "sigmas": [0.2, 0.4, 0.6, 0.8, 1.0],
            },
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": 0,
            "overlap_threshold": 1e-40,
        }

        self.base_molecule_featurizer_test(GMPs, 1)

    def test_features_method_1(self):

        GMPs = {
            "GMPs": {
                "orders": [-1, 0, 1, 2, 3],
                "sigmas": [0.2, 0.4, 0.6, 0.8, 1.0],
            },
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": 1,
            "overlap_threshold": 1e-30,
        }

        self.base_molecule_featurizer_test(GMPs, 1)

    def test_features_method_2(self):

        GMPs = {
            "GMPs": {
                "orders": [-1, 0, 1, 2, 3],
                "sigmas": [0.2, 0.4, 0.6, 0.8, 1.0],
            },
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": 2,
            "overlap_threshold": 1e-30,
        }

        self.base_molecule_featurizer_test(GMPs, 1)

    def test_features_method_3(self):

        GMPs = {
            "GMPs": {
                "orders": [-1, 0, 1, 2, 3],
                "sigmas": [0.2, 0.4, 0.6, 0.8, 1.0],
            },
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": 3,
            "overlap_threshold": 1e-30,
        }

        self.base_molecule_featurizer_test(GMPs, 1)

    def test_features_method_4(self):

        GMPs = {
            "GMPs": {
                "orders": [-1, 0, 1, 2, 3],
                "sigmas": [0.2, 0.4, 0.6, 0.8, 1.0],
            },
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": 4,
            "overlap_threshold": 1e-30,
        }

        self.base_molecule_featurizer_test(GMPs, 1)

    # def test_features_ref_multicore(self):

    #     GMPs = {
    #         "GMPs": {
    #             "orders": [-1, 0, 1, 2, 3],
    #             "sigmas": [0.2, 0.4, 0.6, 0.8, 1.0],
    #         },
    #         "psp_path": "./test_files/NC-SR.gpsp",
    #         "square": False,
    #         "solid_harmonics": True,
    #         "custom_cutoff": -1,
    #         "cutoff": 60,
    #     }

    #     self.base_molecule_featurizer_test(GMPs, 2)

    # def test_features_method_0_multicore(self):

    #     GMPs = {
    #         "GMPs": {
    #             "orders": [-1, 0, 1, 2, 3],
    #             "sigmas": [0.2, 0.4, 0.6, 0.8, 1.0],
    #         },
    #         "psp_path": "./test_files/NC-SR.gpsp",
    #         "square": False,
    #         "solid_harmonics": True,
    #         "custom_cutoff": 0,
    #         "overlap_threshold": 1e-40,
    #     }

    #     self.base_molecule_featurizer_test(GMPs, 2)

    # def test_features_method_1_multicore(self):

    #     GMPs = {
    #         "GMPs": {
    #             "orders": [-1, 0, 1, 2, 3],
    #             "sigmas": [0.2, 0.4, 0.6, 0.8, 1.0],
    #         },
    #         "psp_path": "./test_files/NC-SR.gpsp",
    #         "square": False,
    #         "solid_harmonics": True,
    #         "custom_cutoff": 1,
    #         "overlap_threshold": 1e-30,
    #     }

    #     self.base_molecule_featurizer_test(GMPs, 2)

    # def test_features_method_2_multicore(self):

    #     GMPs = {
    #         "GMPs": {
    #             "orders": [-1, 0, 1, 2, 3],
    #             "sigmas": [0.2, 0.4, 0.6, 0.8, 1.0],
    #         },
    #         "psp_path": "./test_files/NC-SR.gpsp",
    #         "square": False,
    #         "solid_harmonics": True,
    #         "custom_cutoff": 2,
    #         "overlap_threshold": 1e-30,
    #     }

    #     self.base_molecule_featurizer_test(GMPs, 2)

    # def test_features_method_3_multicore(self):

    #     GMPs = {
    #         "GMPs": {
    #             "orders": [-1, 0, 1, 2, 3],
    #             "sigmas": [0.2, 0.4, 0.6, 0.8, 1.0],
    #         },
    #         "psp_path": "./test_files/NC-SR.gpsp",
    #         "square": False,
    #         "solid_harmonics": True,
    #         "custom_cutoff": 3,
    #         "overlap_threshold": 1e-30,
    #     }

    #     self.base_molecule_featurizer_test(GMPs, 2)

    # def test_features_method_4_multicore(self):

    #     GMPs = {
    #         "GMPs": {
    #             "orders": [-1, 0, 1, 2, 3],
    #             "sigmas": [0.2, 0.4, 0.6, 0.8, 1.0],
    #         },
    #         "psp_path": "./test_files/NC-SR.gpsp",
    #         "square": False,
    #         "solid_harmonics": True,
    #         "custom_cutoff": 4,
    #         "overlap_threshold": 1e-30,
    #     }

    #     self.base_molecule_featurizer_test(GMPs, 2)


if __name__ == "__main__":
    unittest.main()
