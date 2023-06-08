# Copyright Toyota Research Institute 2023
"""
Module for defining some of the utility functions needed
"""

import hashlib
import numpy as np
from .constants import ATOM_INDEX_TO_SYMBOL_DICT, ATOM_SYMBOL_TO_INDEX_DICT


def _gen_2Darray_for_ffi(arr, ffi, cdata="double"):
    """
    Function to generate 2D pointer for cffi from 2D numpy array
    Needed to interface the C code
    """
    shape = arr.shape
    arr_p = ffi.new(cdata + " *[%d]" % shape[0])
    for i in range(shape[0]):
        arr_p[i] = ffi.cast(cdata + " *", arr[i].ctypes.data)
    return arr_p


def _gen_2Darray_for_ffi2(arr, ffi, cdata="double"):
    """
    Function to generate 2D pointer for cffi from 2D numpy array
    Needed to interface the C code
    """
    shape = arr.shape
    # dsize = ffi.sizeof(cdata)
    arr_p = ffi.new(cdata + " *[%d]" % shape[0])
    ptr = ffi.cast(cdata + " *", arr.ctypes.data)
    for i in range(shape[0]):
        arr_p[i] = ptr + i * shape[1]
    return arr_p


def get_hash(image, ref_positions):
    """
    Get the hash based on the image object and positions for
    computing GMP features.

    Parameters
    ----------
    image : Dict
        image information in the format of this featrizer,
        containing all information of a system
    ref_positions : numpy.ndarray
        size n*3, Cartesian coordinates of the reference
        points want to be featurized

    Return
    ----------
    complete_hash : str
        hash based on the image and reference points.
    """

    string = ""
    string += str(image["pbc"])

    flattened_cell = image["cell"].flatten()
    for number in np.round(flattened_cell, 4):
        string += "%.4f" % (number + 0)
    for symbol in image["atom_symbols"]:
        string += symbol
    for number in np.round(image["atom_positions"].flatten(), 4):
        string += "%.4f" % (number + 0)
    for number in np.round(image["occupancies"], 4):
        string += "%.4f" % (number + 0)

    md5_1 = hashlib.md5(string.encode("utf-8"))
    hash1 = md5_1.hexdigest()

    pos_string = ""
    for number in np.round(ref_positions.flatten(), 4):
        pos_string += "%.4f" % (number + 0)

    md5_2 = hashlib.md5(pos_string.encode("utf-8"))
    hash2 = md5_2.hexdigest()

    complete_hash = "{}-{}".format(hash1, hash2)
    return complete_hash


def list_symbols_to_indices(list_of_symbols):
    """
    Convert a list of hemical symbols to atomic numbers

    Parameters
    ----------
    list_of_symbols : List
        list of chemical symbols like ["C", "H", ...]

    Return
    ----------
    list_of_indices : List
        list of atomic numbers like [6, 1, ...]
    """

    list_indices = []
    for symbol in list_of_symbols:
        list_indices.append(ATOM_SYMBOL_TO_INDEX_DICT[symbol])
    list_of_indices = np.array(list_indices, dtype=np.intc)
    return list_of_indices


def list_indices_to_symbols(list_of_indices):
    """
    Convert a list of atomic numbers to chemical symbols

    Parameters
    ----------
    list_of_indices : List
        list of atomic numbers like [6, 1, ...]

    Return
    ----------
    list_of_symbols : List
        list of chemical symbols like ["C", "H", ...]
    """

    list_of_symbols = []
    for index in list_of_indices:
        list_of_symbols.append(ATOM_INDEX_TO_SYMBOL_DICT[index])
    return list_of_symbols


# def istarmap(self, func, iterable, chunksize=1):
#     """starmap-version of imap"""
#     self._check_running()
#     if chunksize < 1:
#         raise ValueError("Chunksize must be 1+, not {0:n}".format(chunksize))

#     task_batches = mpp.Pool._get_tasks(func, iterable, chunksize)
#     result = mpp.IMapIterator(self)
#     self._taskqueue.put(
#         (
#             self._guarded_task_generation(result._job, mpp.starmapstar, task_batches),
#             result._set_length,
#         )
#     )
#     return (item for chunk in result for item in chunk)


# def to_iterator(obj_ids):
#     while obj_ids:
#         done, obj_ids = ray.wait(obj_ids)
#         yield ray.get(done[0])


def get_scaled_position(cell, positions):
    """
    Convert atomic positions to scaled atomic positions

    Parameters
    ----------
    cell : numpy.ndarray
        size 3*3, lattice vectors of a unit cell
    positions : numpy.ndarray
        size n*3, positions in Cartesian coordinate

    Return
    ----------
    scaled_positions : numpy.ndarray
        size n*3, corresponding scaled positions
    """

    assert cell.shape == (3, 3)
    scaled_positions = np.linalg.solve(cell.T, np.transpose(positions)).T.round(6)
    return scaled_positions
