import sys
import io
from collections import Counter

# This script requires the 'pymatgen' library.
# We will try to import it and provide a helpful message if it's not found.
try:
    from pymatgen.io.cif import CifParser
    from pymatgen.analysis.local_env import CrystalNN
except ImportError:
    print("Error: The 'pymatgen' library is required to run this script.", file=sys.stderr)
    print("Please install it using the command: pip install pymatgen", file=sys.stderr)
    sys.exit(1)

def find_coordination_polyhedra():
    """
    Parses a CIF string, identifies unique atomic sites, determines their
    coordination polyhedra using pymatgen, and prints the results.
    """
    # The CIF data for the structure provided by the user.
    # Note: The formula from atom counts is ReAl6, not ReAl12.
    cif_data_string = """
data_ReAl6
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
 Al_A      Al     0.3182  0.2158  0.2500
 Al_B      Al     0.0000  0.3662  0.1030
 Al_C      Al     0.1743  0.0000  0.0000
 Re_A      Re     0.0000  0.0445  0.2500
"""

    # 1. Parse the CIF string to create a pymatgen Structure object.
    # The parser uses symmetry operations to build the full unit cell.
    # The `_atom_site_label` is stored as a property on each site.
    parser = CifParser.from_string(cif_data_string)
    structure = parser.get_structures(primitive=False)[0]

    # 2. Initialize the CrystalNN algorithm for coordination analysis.
    cnn = CrystalNN()

    # 3. Get the list of crystallographically unique sites.
    # We only need to analyze one representative site from each Wyckoff position.
    unique_sites = [group[0] for group in structure.equivalent_sites]

    results = []

    # 4. Loop over each unique site to determine its coordination polyhedron.
    for site in unique_sites:
        # Get the label of the central atom (e.g., 'Al_A')
        central_atom_label = site.properties.get('atom_site_label', site.specie.symbol)

        # Find the nearest neighbors using CrystalNN
        nn_info = cnn.get_nn_info(structure, site.index)

        # Count the neighbors by their element type
        neighbor_elements = [nn['site'].specie.symbol for nn in nn_info]
        neighbor_counts = Counter(neighbor_elements)

        # 5. Build the chemical formula string for the coordination polyhedron.
        formula_parts = []
        # Sort by element symbol for a consistent formula order (Al, then Re)
        for element in sorted(neighbor_counts.keys()):
            count = neighbor_counts[element]
            # Append element and its count. Omit count if it's 1.
            formula_parts.append(f"{element}{count if count > 1 else ''}")
        poly_formula = "".join(formula_parts)

        # Format the result string for this site
        results.append(f"{central_atom_label}, {poly_formula}")

    # 6. Sort the final results alphabetically by the central atom's label
    # This ensures a consistent and predictable output order.
    results.sort()

    # 7. Print the final answer in the required format.
    print("; ".join(results))

if __name__ == "__main__":
    find_coordination_polyhedra()