import pymatgen.core as mg
from pymatgen.analysis.local_env import VoronoiNN
from collections import Counter
import io

def find_coordination_polyhedra():
    """
    Parses a CIF string, finds the coordination environment for each unique site,
    and prints the results in the specified format.
    """
    cif_data = """
_cell_length_a    7.609100
_cell_length_b    6.611700
_cell_length_c    9.023000
_cell_angle_alpha 90.000000
_cell_angle_beta  90.000000
_cell_angle_gamma 90.000000
_symmetry_space_group_name_H-M Cmcm

loop_
 _atom_site_label
 _atom_site_type_symbol
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 Al_A      Al  0.3182  0.2158  0.2500
 Al_B      Al  0.0000  0.3662  0.1030
 Al_C      Al  0.1743  0.0000  0.0000
 Re_A      Re  0.0000  0.0445  0.2500
"""
    
    # 1. Parse the structure from the CIF string
    structure = mg.Structure.from_file(io.StringIO(cif_data), fmt="cif")

    # 2. Initialize the nearest neighbor finder using the Voronoi method
    nn_finder = VoronoiNN()

    # 3. Store results to avoid duplicates from symmetry-equivalent sites
    results = set()

    # 4. Iterate through each unique site in the structure
    for i in range(len(structure)):
        central_atom_symbol = structure[i].specie.symbol

        # Get the list of nearest neighbors
        neighbors_info = nn_finder.get_nn_info(structure, i)
        
        # Count the number of neighbors of each element type
        neighbor_elements = [neighbor['site'].specie.symbol for neighbor in neighbors_info]
        counts = Counter(neighbor_elements)
        
        # Build the formula string for the polyhedron
        # The elements in the formula are sorted alphabetically for consistency
        formula_parts = []
        for elem in sorted(counts.keys()):
            count = counts[elem]
            # Omit the number '1' for single atoms
            if count == 1:
                formula_parts.append(elem)
            else:
                formula_parts.append(f"{elem}{count}")
        
        formula = "".join(formula_parts)
        
        # Add the 'central_atom, formula' string to the results set
        results.add(f"{central_atom_symbol}, {formula}")

    # Sort the results alphabetically and print them in the required format
    print("; ".join(sorted(list(results))))

if __name__ == "__main__":
    find_coordination_polyhedra()