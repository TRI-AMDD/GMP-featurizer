import multiprocessing.pool as mpp
import ray
import hashlib

import numpy as np

from .constants import ATOM_INDEX_TO_SYMBOL_DICT, ATOM_SYMBOL_TO_INDEX_DICT

# from ase.io.trajectory import Trajectory


def _gen_2Darray_for_ffi(arr, ffi, cdata="double"):
    # Function to generate 2D pointer for cffi
    shape = arr.shape
    arr_p = ffi.new(cdata + " *[%d]" % shape[0])
    for i in range(shape[0]):
        arr_p[i] = ffi.cast(cdata + " *", arr[i].ctypes.data)
    return arr_p


def _gen_2Darray_for_ffi2(arr, ffi, cdata="double"):
    # Function to generate 2D pointer for cffi
    shape = arr.shape
    dsize = ffi.sizeof(cdata)
    arr_p = ffi.new(cdata + " *[%d]" % shape[0])
    ptr = ffi.cast(cdata + " *", arr.ctypes.data)
    for i in range(shape[0]):
        arr_p[i] = ptr + i * shape[1]
    return arr_p


# def get_hash(image, ref_positions):

#     string = ""
#     string += str(image.pbc)
#     try:
#         flattened_cell = image.cell.array.flatten()
#     except AttributeError:  # older ASE
#         flattened_cell = image.cell.flatten()
#     for number in flattened_cell:
#         string += "%.15f" % number
#     for number in image.get_atomic_numbers():
#         string += "%3d" % number
#     for number in image.get_positions(wrap=True).flatten():
#         string += "%.15f" % number

#     md5_1 = hashlib.md5(string.encode("utf-8"))
#     hash1 = md5_1.hexdigest()

#     pos_string = ""
#     for number in ref_positions.flatten():
#         pos_string += "%.15f" % number

#     md5_2 = hashlib.md5(pos_string.encode("utf-8"))
#     hash2 = md5_2.hexdigest()

#     return "{}-{}".format(hash1, hash2)


def get_hash(image, ref_positions):

    string = ""
    string += str(image["pbc"])

    flattened_cell = image["cell"].flatten()
    for number in flattened_cell:
        string += "%.15f" % number
    for symbol in image["atom_symbols"]:
        string += symbol
    for number in image["atom_positions"].flatten():
        string += "%.15f" % number
    for number in image["occupancies"]:
        string += "%.15f" % number

    md5_1 = hashlib.md5(string.encode("utf-8"))
    hash1 = md5_1.hexdigest()

    pos_string = ""
    for number in ref_positions.flatten():
        pos_string += "%.15f" % number

    md5_2 = hashlib.md5(pos_string.encode("utf-8"))
    hash2 = md5_2.hexdigest()

    return "{}-{}".format(hash1, hash2)


# def validate_image(image, ref_positions):
#     for number in image.get_scaled_positions(wrap=True).flatten():
#         if number > 1.0 or number < 0.0:
#             raise ValueError(
#                 "****ERROR: scaled position not strictly between [0, 1]"
#                 "Please check atom position and system cell size are set up correctly"
#             )
#     return


def list_symbols_to_indices(list_of_symbols):
    list_indices = []
    for symbol in list_of_symbols:
        list_indices.append(ATOM_SYMBOL_TO_INDEX_DICT[symbol])
    return np.array(list_indices, dtype=np.intc)


def list_indices_to_symbols(list_of_indices):
    list_symbols = []
    for index in list_of_indices:
        list_symbols.append(ATOM_INDEX_TO_SYMBOL_DICT[index])
    return list_symbols


def istarmap(self, func, iterable, chunksize=1):
    """starmap-version of imap"""
    self._check_running()
    if chunksize < 1:
        raise ValueError("Chunksize must be 1+, not {0:n}".format(chunksize))

    task_batches = mpp.Pool._get_tasks(func, iterable, chunksize)
    result = mpp.IMapIterator(self)
    self._taskqueue.put(
        (
            self._guarded_task_generation(result._job, mpp.starmapstar, task_batches),
            result._set_length,
        )
    )
    return (item for chunk in result for item in chunk)


def to_iterator(obj_ids):
    while obj_ids:
        done, obj_ids = ray.wait(obj_ids)
        yield ray.get(done[0])


def get_scaled_position(cell, positions):
    assert cell.shape == (3, 3)
    return np.linalg.solve(cell.T, np.transpose(positions)).T.round(15)
