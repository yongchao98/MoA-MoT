import numpy as np
from collections import Counter
import math

def find_coordination_polyhedra():
    """
    This script analyzes the provided CIF data for ReAl12 to determine
    the coordination polyhedra for each unique atomic site.
    """
    # 1. Parse Crystal Data from CIF
    cell_lengths = {'a': 7.6091, 'b': 6.6117, 'c': 9.0230}
    cell_angles = {'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0}

    unique_atoms = [
        {'label': 'Al_A', 'type': 'Al', 'coords': np.array([0.3182, 0.2158, 0.2500])},
        {'label': 'Al_B', 'type': 'Al', 'coords': np.array([0.0000, 0.3662, 0.1030])},
        {'label': 'Al_C', 'type': 'Al', 'coords': np.array([0.1743, 0.0000, 0.0000])},
        {'label': 'Re_A', 'type': 'Re', 'coords': np.array([0.0000, 0.0445, 0.2500])},
    ]

    sym_ops_str = [
        '-x+1/2,-y+1/2,-z', ' -x+1/2,-y+1/2,z+1/2', ' -x+1/2,y+1/2,-z+1/2',
        '-x+1/2,y+1/2,z', ' -x,-y,-z', ' -x,-y,z+1/2', ' -x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', ' x+1/2,-y+1/2,z+1/2', ' x+1/2,y+1/2,-z+1/2',
        'x+1/2,y+1/2,z', ' x,-y,-z', ' x,-y,z+1/2', ' x,y,-z+1/2', 'x,y,z'
    ]

    # Helper function to parse symmetry operation strings
    def parse_sym_op(op_str):
        return lambda x, y, z: tuple(eval(s) for s in op_str.split(','))

    sym_ops = [parse_sym_op(s) for s in sym_ops_str]

    # 2. Generate Full Unit Cell
    unit_cell_atoms = []
    seen_coords = set()
    for atom in unique_atoms:
        for op in sym_ops:
            x, y, z = atom['coords']
            new_x, new_y, new_z = op(x, y, z)
            # Normalize coordinates to be within [0, 1)
            new_coords = np.array([
                new_x - math.floor(new_x),
                new_y - math.floor(new_y),
                new_z - math.floor(new_z)
            ])
            # Round to avoid floating point duplicates
            coord_tuple = tuple(np.round(new_coords, 4))
            if coord_tuple not in seen_coords:
                unit_cell_atoms.append({'type': atom['type'], 'coords': new_coords})
                seen_coords.add(coord_tuple)

    # Build fractional-to-Cartesian transformation matrix
    a, b, c = cell_lengths['a'], cell_lengths['b'], cell_lengths['c']
    # For orthorhombic system, the matrix is diagonal
    frac_to_cart_matrix = np.array([[a, 0, 0], [0, b, 0], [0, 0, c]])

    # 3. Construct a Supercell
    supercell_atoms = []
    for atom in unit_cell_atoms:
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    supercell_atoms.append({
                        'type': atom['type'],
                        'coords': atom['coords'] + np.array([i, j, k])
                    })

    polyhedra_results = set()

    # 4-6. Iterate through unique atoms to find coordination
    for central_atom_info in unique_atoms:
        central_type = central_atom_info['type']
        central_frac_coords = central_atom_info['coords']

        # Calculate distances to all atoms in the supercell
        distances = []
        for neighbor_atom in supercell_atoms:
            vec_frac = central_frac_coords - neighbor_atom['coords']
            vec_cart = np.dot(frac_to_cart_matrix, vec_frac)
            dist = np.linalg.norm(vec_cart)
            
            if dist > 1e-4:  # Exclude the atom itself
                distances.append({'dist': dist, 'type': neighbor_atom['type']})
        
        distances.sort(key=lambda x: x['dist'])

        if not distances:
            continue
            
        # Identify the first coordination shell by finding the largest gap in distances
        unique_dists = sorted(list(set(d['dist'] for d in distances)))
        cutoff = unique_dists[0] * 1.2 # Default cutoff if only one shell found
        if len(unique_dists) > 1:
            gaps = [unique_dists[i] - unique_dists[i-1] for i in range(1, len(unique_dists))]
            max_gap_index = np.argmax(gaps)
            # Set cutoff to be the distance just before the largest gap
            cutoff = unique_dists[max_gap_index] + 1e-3

        first_shell_atoms = [d['type'] for d in distances if d['dist'] <= cutoff]

        # Determine the formula of the polyhedron
        counts = Counter(first_shell_atoms)
        formula_parts = []
        for elem in sorted(counts.keys()): # Sort for consistent formula (Al, Re)
            count = counts[elem]
            formula_parts.append(f"{elem}{count if count > 1 else ''}")
        formula = "".join(formula_parts)
        
        polyhedra_results.add(f"{central_type}, {formula}")

    # 7. Format and print the final output
    final_output = "; ".join(sorted(list(polyhedra_results)))
    print(final_output)

if __name__ == '__main__':
    find_coordination_polyhedra()