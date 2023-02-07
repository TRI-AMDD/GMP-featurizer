import unittest
import numpy as np
from GMPFeaturizer import GMPFeaturizer
import pickle


class FeatruizerTest(unittest.TestCase):
    def test_features_molecules_ref(self):
        # ensuring the reference method give expected outputs
        with open("./test_files/molecules.p", "rb") as f:
            images = pickle.load(f)

        with open("./test_files/molecules_ref_features.p", "rb") as f:
            reference_features = pickle.load(f)

        # GMP hyperparameters
        GMP_order = 3
        nsigmas = 5
        width = 1.0

        # setup featurizer/calcualtor
        sigmas = np.round(np.linspace(0.0, width, nsigmas + 1, endpoint=True), 4)[1:]
        GMPs = {
            "GMPs": {"orders": [-1] + list(range(GMP_order + 1)), "sigmas": sigmas},
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": -1,
            "cutoff": 60,
        }

        featurizer_noderiv_ref = GMPFeaturizer(
            GMPs=GMPs, calc_derivatives=False, calc_occ_derivatives=False
        )
        features_noderiv_ref = featurizer_noderiv_ref.prepare_features(images)

        featurizer_fp_deriv_ref = GMPFeaturizer(
            GMPs=GMPs, calc_derivatives=True, calc_occ_derivatives=False
        )
        features_fp_deriv_ref = featurizer_fp_deriv_ref.prepare_features(images)

        featurizer_occ_deriv_ref = GMPFeaturizer(
            GMPs=GMPs, calc_derivatives=False, calc_occ_derivatives=True
        )
        features_occ_deriv_ref = featurizer_occ_deriv_ref.prepare_features(images)

        featurizer_fp_occ_deriv_ref = GMPFeaturizer(
            GMPs=GMPs, calc_derivatives=True, calc_occ_derivatives=True
        )
        features_fp_occ_deriv_ref = featurizer_fp_occ_deriv_ref.prepare_features(images)

        assert np.all(
            [
                np.allclose(
                    features_noderiv_ref[i]["features"],
                    reference_features[i]["features"],
                    rtol=1e-12,
                    atol=1e-18,
                )
                for i in range(len(images))
            ]
        )
        assert np.all(
            [
                np.allclose(
                    features_fp_deriv_ref[i]["features"],
                    reference_features[i]["features"],
                    rtol=1e-12,
                    atol=1e-18,
                )
                for i in range(len(images))
            ]
        )
        assert np.all(
            [
                np.allclose(
                    features_occ_deriv_ref[i]["features"],
                    reference_features[i]["features"],
                    rtol=1e-12,
                    atol=1e-18,
                )
                for i in range(len(images))
            ]
        )
        assert np.all(
            [
                np.allclose(
                    features_fp_occ_deriv_ref[i]["features"],
                    reference_features[i]["features"],
                    rtol=1e-12,
                    atol=1e-18,
                )
                for i in range(len(images))
            ]
        )

    def test_features_molecules_other(self):
        # ensuring all the methods give expected outputs
        with open("./test_files/molecules.p", "rb") as f:
            images = pickle.load(f)

        with open("./test_files/molecules_ref_features.p", "rb") as f:
            reference_features = pickle.load(f)

        # GMP hyperparameters
        GMP_order = 3
        nsigmas = 5
        width = 1.0

        # setup featurizer/calcualtor
        sigmas = np.round(np.linspace(0.0, width, nsigmas + 1, endpoint=True), 4)[1:]
        for setting in [0, 1, 2, 3, 4]:
            GMPs = {
                "GMPs": {"orders": [-1] + list(range(GMP_order + 1)), "sigmas": sigmas},
                "psp_path": "./test_files/NC-SR.gpsp",
                "square": False,
                "solid_harmonics": True,
                "custom_cutoff": setting,
                "overlap_threshold": 1e-30,
            }

            featurizer_noderiv = GMPFeaturizer(
                GMPs=GMPs, calc_derivatives=False, calc_occ_derivatives=False
            )
            features_noderiv = featurizer_noderiv.prepare_features(images)

            featurizer_fp_deriv = GMPFeaturizer(
                GMPs=GMPs, calc_derivatives=True, calc_occ_derivatives=False
            )
            features_fp_deriv = featurizer_fp_deriv.prepare_features(images)

            featurizer_occ_deriv = GMPFeaturizer(
                GMPs=GMPs, calc_derivatives=False, calc_occ_derivatives=True
            )
            features_occ_deriv = featurizer_occ_deriv.prepare_features(images)

            featurizer_fp_occ_deriv = GMPFeaturizer(
                GMPs=GMPs, calc_derivatives=True, calc_occ_derivatives=True
            )
            features_fp_occ_deriv = featurizer_fp_occ_deriv.prepare_features(images)

            assert np.all(
                [
                    np.allclose(
                        features_noderiv[i]["features"],
                        reference_features[i]["features"],
                        rtol=1e-12,
                        atol=1e-18,
                    )
                    for i in range(len(images))
                ]
            )
            assert np.all(
                [
                    np.allclose(
                        features_fp_deriv[i]["features"],
                        reference_features[i]["features"],
                        rtol=1e-12,
                        atol=1e-18,
                    )
                    for i in range(len(images))
                ]
            )
            assert np.all(
                [
                    np.allclose(
                        features_occ_deriv[i]["features"],
                        reference_features[i]["features"],
                        rtol=1e-12,
                        atol=1e-18,
                    )
                    for i in range(len(images))
                ]
            )
            assert np.all(
                [
                    np.allclose(
                        features_fp_occ_deriv[i]["features"],
                        reference_features[i]["features"],
                        rtol=1e-12,
                        atol=1e-18,
                    )
                    for i in range(len(images))
                ]
            )


if __name__ == "__main__":
    unittest.main()
