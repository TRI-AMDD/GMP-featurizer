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
    """
    shape = arr.shape
    arr_p = ffi.new(cdata + " *[%d]" % shape[0])
    for i in range(shape[0]):
        arr_p[i] = ffi.cast(cdata + " *", arr[i].ctypes.data)
    return arr_p


def _gen_2Darray_for_ffi2(arr, ffi, cdata="double"):
    """
    Function to generate 2D pointer for cffi from 2D numpy array
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
    """

    string = ""
    string += str(image["pbc"])

    flattened_cell = image["cell"].flatten()
    for number in flattened_cell:
        string += "%.5f" % number
    for symbol in image["atom_symbols"]:
        string += symbol
    for number in image["atom_positions"].flatten():
        string += "%.5f" % number
    for number in image["occupancies"]:
        string += "%.5f" % number

    md5_1 = hashlib.md5(string.encode("utf-8"))
    hash1 = md5_1.hexdigest()

    pos_string = ""
    for number in ref_positions.flatten():
        pos_string += "%.5f" % number

    md5_2 = hashlib.md5(pos_string.encode("utf-8"))
    hash2 = md5_2.hexdigest()

    return "{}-{}".format(hash1, hash2)


def list_symbols_to_indices(list_of_symbols):
    """
    Chemical symbols to atomic number
    """
    list_indices = []
    for symbol in list_of_symbols:
        list_indices.append(ATOM_SYMBOL_TO_INDEX_DICT[symbol])
    return np.array(list_indices, dtype=np.intc)


def list_indices_to_symbols(list_of_indices):
    """
    Atomic number to chemical symbols
    """
    list_symbols = []
    for index in list_of_indices:
        list_symbols.append(ATOM_INDEX_TO_SYMBOL_DICT[index])
    return list_symbols


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
    atomic positions to scaled atomic positions
    """
    assert cell.shape == (3, 3)
    return np.linalg.solve(cell.T, np.transpose(positions)).T.round(6)
