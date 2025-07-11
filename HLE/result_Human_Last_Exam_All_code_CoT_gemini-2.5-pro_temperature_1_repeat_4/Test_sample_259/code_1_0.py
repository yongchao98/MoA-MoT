import pymatgen.core as mg
from pymatgen.analysis.local_env import CrystalNN
import collections

def find_coordination_polyhedra():
    """
    Parses a CIF string, identifies the coordination environment for each unique
    atom site, and prints the results in the specified format.
    """
    # The CIF data for the structure is provided here as a multiline string.
    # Note: The atomic coordinates and stoichiometry correspond to ReAl6,
    # and the analysis is performed on this structure.
    cif_data = """
data_ReAl6
_cell_length_a    7.609100
_cell_length_b    6.611700
_cell_length_c    9.023000
_cell_angle_alpha 90.000000
_cell_angle_beta  90.000000
_cell_angle_gamma 90.000000
_symmetry_space_group_name_H-M 'Cmcm'

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

    try:
        # Create a Pymatgen Structure object from the CIF data string.
        structure = mg.Structure.from_str(cif_data, fmt="cif")

        # Initialize the CrystalNN algorithm to determine coordination environments.
        cnn = CrystalNN()

        # A list to hold the formatted strings for each polyhedron.
        results = []

        # Iterate through each unique site in the asymmetric unit of the structure.
        for i, site in enumerate(structure.sites):
            # Get the symbol of the central atom (e.g., 'Al', 'Re').
            central_atom_symbol = site.specie.symbol

            # Use CrystalNN to find the nearest neighbors.
            neighbors = cnn.get_cn_neighbors(structure, n=i)

            # Count the neighbors by element type.
            neighbor_counts = collections.Counter([s.specie.symbol for s in neighbors])
            
            # Create the polyhedron formula string, e.g., Al10Re2
            # The elements are sorted alphabetically for consistent representation.
            poly_formula = "".join(
                f"{elem}{neighbor_counts[elem]}"
                for elem in sorted(neighbor_counts)
            )

            # Format the output string as "central atom, formula".
            result_string = f"{central_atom_symbol}, {poly_formula}"
            results.append(result_string)

        # Join all the individual results with a semicolon and a space.
        print("; ".join(results))

    except ImportError:
        print("Error: This script requires the 'pymatgen' library.")
        print("Please install it using: pip install pymatgen")
    except Exception as e:
        print(f"An error occurred: {e}")

# Execute the function to find and print the coordination polyhedra.
find_coordination_polyhedra()