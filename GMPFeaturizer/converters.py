import numpy as np


class ASEAtomsConverter():
    def __init__(self):
        pass
    
    def convert(self, atoms_list):
        results = []
        for atoms in atoms_list:
            temp = {}
            temp["cell"]=atoms.get_cell()[:]
            temp["pbc"]=np.copy(atoms.get_pbc()).astype(np.intc)
            temp["atom_positions"]=atoms.get_positions(wrap=True)
            temp["atom_symbols"]=atoms.get_chemical_symbols()
            temp["occupancies"]=np.array([1.0 for _ in range(len(atoms.get_chemical_symbols()))])
            results.append(temp)
        return results


class PymatgenStructureConverter():
    def __init__(self):
        pass
    
    def convert(self, structure_list):
        results = []
        for structure in structures_list:
            temp = {}
            temp["cell"] = np.array(structure.lattice.matrix)
            # ref_positions = np.array([site.coords for site in structure.sites])
            atom_positions = []
            atom_symbols = []
            occupancies = []
            for site in structure.sites:
                for el, occ in site.species.fractional_composition.get_el_amt_dict().items():
                    atom_positions.append(site.coords)
                    atom_symbols.append(el)
                    occupancies.append(occ)
            temp["atom_positions"] = np.array(atom_positions)
            temp["atom_symbols"] = atom_symbols
            temp["occupancies"] = np.array(occupancies)
            temp["pbc"] = np.copy(np.array(structure.lattice.pbc)).astype(np.intc)
            results.append(temp)
        return results

# ref_positions = image.get_positions()
# symbols = np.array(image.get_chemical_symbols())
# atom_num = len(symbols)
# # atom_indices = list_symbols_to_indices(symbols)
# cell = image.cell
# scaled_ref_positions = cell.scaled_positions(ref_positions)
# scaled_ref_positions = np.array(
#     [np.array(v, dtype="float64") for v in scaled_ref_positions],
#     dtype="float64",
# )
# # atom_indices_p = ffi.cast("int *", atom_indices.ctypes.data)

# cart = np.copy(image.get_positions(wrap=True), order="C")
# scale = np.copy(image.get_scaled_positions(wrap=True), order="C")
# cell = np.copy(image.cell, order="C")
# pbc = np.copy(image.get_pbc()).astype(np.intc)