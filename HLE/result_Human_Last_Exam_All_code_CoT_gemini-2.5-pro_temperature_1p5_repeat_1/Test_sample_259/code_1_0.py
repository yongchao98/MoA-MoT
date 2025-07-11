import sys
try:
    import pymatgen.core as mg
    from pymatgen.analysis.local_env import VoronoiNN
except ImportError:
    # If pymatgen is not installed, print an error message to stderr and exit.
    print("Error: The 'pymatgen' library is required to run this script.", file=sys.stderr)
    print("Please install it using: pip install pymatgen", file=sys.stderr)
    sys.exit(1)

# The CIF data for ReAl12 provided in the problem description.
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

def find_coordination_polyhedra():
    """
    Parses CIF data, finds coordination polyhedra for each unique site,
    and prints the result in the specified format.
    """
    # 1. Parse the CIF data into a Pymatgen Structure object.
    structure = mg.Structure.from_str(cif_data, fmt="cif")

    # 2. Initialize the Voronoi-based nearest neighbor finding algorithm.
    nn_finder = VoronoiNN()

    # 3. Iterate through each unique site to find its coordination.
    # The structure object from the CIF contains only unique crystallographic sites.
    output_parts = []
    # Sites are read in the order they appear in the CIF file.
    # We will process Al_A, Al_B, Al_C, then Re_A.
    for i, site in enumerate(structure):
        # Get the chemical symbol of the central atom.
        central_atom_symbol = site.specie.symbol

        # Get a dictionary of coordinating atoms and their counts, e.g., {'Al': 12, 'Re': 2}.
        coord_dict = nn_finder.get_cn_dict(structure, i)

        # 4. Construct the chemical formula for the coordination polyhedron.
        formula_parts = []
        # Sort by element symbol for a consistent and readable formula.
        for element_symbol in sorted(coord_dict.keys()):
            count = coord_dict[element_symbol]
            # Append element and its count. If count is 1, omit the number.
            formula_parts.append(f"{element_symbol}{count if count > 1 else ''}")
        
        polyhedron_formula = "".join(formula_parts)

        # Combine the central atom and the polyhedron formula.
        output_parts.append(f"{central_atom_symbol}, {polyhedron_formula}")

    # 5. Print the final results in the specified format.
    print("; ".join(output_parts))

if __name__ == "__main__":
    find_coordination_polyhedra()