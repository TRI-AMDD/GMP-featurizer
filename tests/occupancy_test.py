import unittest
import numpy as np
from GMPFeaturizer import GMPFeaturizer
from GMPFeaturizer import ASEAtomsConverter
import pickle


class FeatruizerOccTest(unittest.TestCase):
    def base_occupany_derivative_test(self, GMPs, atom_index):
        # ensuring the reference method give expected outputs
        with open("./test_files/molecules.p", "rb") as f:
            images = pickle.load(f)

        converter = ASEAtomsConverter()

        print("=========================")
        print("test atom {}".format(atom_index))
        print("test 1 (both derivatives calcualted)")
        featurizer = GMPFeaturizer(
            GMPs=GMPs, calc_derivatives=True, calc_occ_derivatives=True, verbose=False
        )
        print("docc\tmax difference\t\ttolerance\tpassed")

        for docc, tol in [
            (0.002, 1e-4),
            (0.005, 3e-4),
            (0.01, 1e-3),
            (0.02, 3e-3),
            (0.05, 1e-2),
            (0.1, 3e-2),
            (0.2, 1e-1),
        ]:
            converted_images = converter.convert(images)
            converted_images2 = converter.convert(images)
            for image in converted_images2:
                image["occupancies"][atom_index] = 1 - docc
            features1 = featurizer.prepare_features(converted_images, cores=1)
            features2 = featurizer.prepare_features(converted_images2, cores=1)

            max_error_list = []
            for i in range(len(images)):
                ref = features1[i]["features"]
                pred = (
                    features2[i]["feature_occ_primes"][:, atom_index].reshape(
                        features2[i]["features"].shape
                    )
                    * docc
                    + features2[i]["features"]
                )
                max_error_list.append(np.max(np.abs(ref - pred)))
            max_error = np.max(max_error_list)
            print("{}\t{}\t{}\t\t{}".format(docc, max_error, tol, max_error < tol))
            assert max_error < tol

        print("------")
        print("test 2 (only occ derivatives calcualted)")
        featurizer = GMPFeaturizer(
            GMPs=GMPs, calc_derivatives=False, calc_occ_derivatives=True, verbose=False
        )
        print("docc\tmax difference\t\ttolerance\tpassed")

        for docc, tol in [
            (0.002, 1e-4),
            (0.005, 3e-4),
            (0.01, 1e-3),
            (0.02, 3e-3),
            (0.05, 1e-2),
            (0.1, 3e-2),
            (0.2, 1e-1),
        ]:
            converted_images = converter.convert(images)
            converted_images2 = converter.convert(images)
            for image in converted_images2:
                image["occupancies"][atom_index] = 1 - docc
            features1 = featurizer.prepare_features(converted_images, cores=1)
            features2 = featurizer.prepare_features(converted_images2, cores=1)

            max_error_list = []
            for i in range(len(images)):
                ref = features1[i]["features"]
                pred = (
                    features2[i]["feature_occ_primes"][:, atom_index].reshape(
                        features2[i]["features"].shape
                    )
                    * docc
                    + features2[i]["features"]
                )
                max_error_list.append(np.max(np.abs(ref - pred)))
            max_error = np.max(max_error_list)
            print("{}\t{}\t{}\t\t{}".format(docc, max_error, tol, max_error < tol))
            assert max_error < tol

    def test_occupancy_drivative_ref(self):

        GMPs = {
            "GMPs": {"orders": [-1, 0, 1, 2], "sigmas": [0.1, 0.2, 0.3],},
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": -1,
            "cutoff": 60,
        }

        self.base_occupany_derivative_test(GMPs, 0)
        self.base_occupany_derivative_test(GMPs, 1)
        self.base_occupany_derivative_test(GMPs, 2)
        self.base_occupany_derivative_test(GMPs, 3)
        self.base_occupany_derivative_test(GMPs, 4)

    def test_occupancy_drivative_method_0(self):

        GMPs = {
            "GMPs": {"orders": [-1, 0, 1, 2], "sigmas": [0.1, 0.2, 0.3],},
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": 0,
            "overlap_threshold": 1e-40,
        }

        self.base_occupany_derivative_test(GMPs, 0)
        self.base_occupany_derivative_test(GMPs, 1)
        self.base_occupany_derivative_test(GMPs, 2)
        self.base_occupany_derivative_test(GMPs, 3)
        self.base_occupany_derivative_test(GMPs, 4)

    def test_occupancy_drivative_method_1(self):

        GMPs = {
            "GMPs": {"orders": [-1, 0, 1, 2], "sigmas": [0.1, 0.2, 0.3],},
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": 1,
            "overlap_threshold": 1e-30,
        }

        self.base_occupany_derivative_test(GMPs, 0)
        self.base_occupany_derivative_test(GMPs, 1)
        self.base_occupany_derivative_test(GMPs, 2)
        self.base_occupany_derivative_test(GMPs, 3)
        self.base_occupany_derivative_test(GMPs, 4)

    def test_occupancy_drivative_method_2(self):

        GMPs = {
            "GMPs": {"orders": [-1, 0, 1, 2], "sigmas": [0.1, 0.2, 0.3],},
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": 2,
            "overlap_threshold": 1e-30,
        }

        self.base_occupany_derivative_test(GMPs, 0)
        self.base_occupany_derivative_test(GMPs, 1)
        self.base_occupany_derivative_test(GMPs, 2)
        self.base_occupany_derivative_test(GMPs, 3)
        self.base_occupany_derivative_test(GMPs, 4)

    def test_occupancy_drivative_method_3(self):

        GMPs = {
            "GMPs": {"orders": [-1, 0, 1, 2], "sigmas": [0.1, 0.2, 0.3],},
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": 3,
            "overlap_threshold": 1e-30,
        }

        self.base_occupany_derivative_test(GMPs, 0)
        self.base_occupany_derivative_test(GMPs, 1)
        self.base_occupany_derivative_test(GMPs, 2)
        self.base_occupany_derivative_test(GMPs, 3)
        self.base_occupany_derivative_test(GMPs, 4)

    def test_occupancy_drivative_method_4(self):

        GMPs = {
            "GMPs": {"orders": [-1, 0, 1, 2], "sigmas": [0.1, 0.2, 0.3],},
            "psp_path": "./test_files/NC-SR.gpsp",
            "square": False,
            "solid_harmonics": True,
            "custom_cutoff": 4,
            "overlap_threshold": 1e-30,
        }

        self.base_occupany_derivative_test(GMPs, 0)
        self.base_occupany_derivative_test(GMPs, 1)
        self.base_occupany_derivative_test(GMPs, 2)
        self.base_occupany_derivative_test(GMPs, 3)
        self.base_occupany_derivative_test(GMPs, 4)


if __name__ == "__main__":
    unittest.main()
