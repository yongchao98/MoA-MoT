import io
from collections import Counter

# To run this code, you need to install the pymatgen library:
# pip install pymatgen

try:
    from pymatgen.core import Structure
    from pymatgen.analysis.local_env import VoronoiNN
except ImportError:
    print("Pymatgen library not found. Please install it using 'pip install pymatgen'")
    exit()

# The CIF data for ReAl12 provided by the user
cif_data = """
data_ReAl12
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

def get_coordination_polyhedra(cif_string):
    """
    Parses a CIF string, finds coordination polyhedra for each unique site,
    and returns the result in the specified format.
    """
    # Use io.StringIO to treat the string as a file for pymatgen
    structure = Structure.from_str(cif_string, fmt="cif")

    # Initialize the Voronoi nearest neighbor finder.
    # A small tolerance can help ignore artifacts from tiny Voronoi facets.
    finder = VoronoiNN(tol=0.2)

    # The input CIF contains the asymmetric unit, so each site in the
    # initial structure object is crystallographically unique.
    unique_sites_indices = range(len(structure))

    results = []
    # Loop through each unique site to determine its coordination environment
    for i in unique_sites_indices:
        central_site = structure[i]
        central_atom_symbol = central_site.specie.symbol
        
        # Get a list of nearest neighbors
        neighbors_info = finder.get_nn_info(structure, i)
        
        # Extract the element symbols of the neighbors
        neighbor_symbols = [info['site'].specie.symbol for info in neighbors_info]
        
        # Count the number of neighbors of each element type
        counts = Counter(neighbor_symbols)
        
        # Build the chemical formula string for the polyhedron
        formula_parts = []
        # Sort by element symbol for a consistent formula (e.g., Al before Re)
        for symbol in sorted(counts.keys()):
            count = counts[symbol]
            formula_parts.append(f"{symbol}{count if count > 1 else ''}")
        
        poly_formula = "".join(formula_parts)
        
        results.append(f"{central_atom_symbol}, {poly_formula}")

    # Print the final result string
    print("; ".join(results))

# Run the analysis and print the result
get_coordination_polyhedra(cif_data)