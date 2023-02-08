import numpy as np
from abc import ABC, abstractmethod


class ImageObjectConverter:
    """
    Base class for image object converters
    """

    def __init__(self):
        pass

    @abstractmethod
    def convert(self, object_list):
        pass


class ASEAtomsConverter(ImageObjectConverter):
    """
    Converter for converting list of ase atoms to a list
    of information that can be read by the featurizer
    """

    def convert(self, object_list):
        results = []
        for atoms in object_list:
            temp = {}
            temp["cell"] = atoms.get_cell()[:]
            temp["pbc"] = np.copy(atoms.get_pbc()).astype(np.intc)
            temp["atom_positions"] = atoms.get_positions(wrap=True)
            temp["atom_symbols"] = atoms.get_chemical_symbols()
            temp["occupancies"] = np.array(
                [1.0 for _ in range(len(atoms.get_chemical_symbols()))]
            )
            results.append(temp)
        return results


class PymatgenStructureConverter(ImageObjectConverter):
    """
    Converter for converting list of pymatgen structures to a list
    of information that can be read by the featurizer
    """

    def convert(self, object_list):
        results = []
        for structure in object_list:
            temp = {}
            temp["cell"] = np.array(structure.lattice.matrix)
            # ref_positions = np.array([site.coords for site in structure.sites])
            atom_positions = []
            atom_symbols = []
            occupancies = []
            for site in structure.sites:
                for (
                    el,
                    occ,
                ) in site.species.fractional_composition.get_el_amt_dict().items():
                    atom_positions.append(site.coords)
                    atom_symbols.append(el)
                    occupancies.append(occ)
            temp["atom_positions"] = np.array(atom_positions)
            temp["atom_symbols"] = atom_symbols
            temp["occupancies"] = np.array(occupancies)
            temp["pbc"] = np.copy(np.array(structure.lattice.pbc)).astype(np.intc)
            results.append(temp)
        return results
