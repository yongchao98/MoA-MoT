import numpy as np
from collections import Counter

def find_coordination_polyhedra():
    """
    Analyzes a crystal structure defined in the script to find the
    coordination polyhedra for each unique atomic site.
    """
    # Part 1: Define the crystal structure from the CIF file
    cell_lengths = np.array([7.609100, 6.611700, 9.023000])
    
    # Symmetry operations for space group Cmcm (No. 63)
    sym_ops_str = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2',
        '-x+1/2,y+1/2,z', '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2',
        'x+1/2,y+1/2,z', 'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]

    # Atoms in the asymmetric unit
    au_atoms = [
        {'label': 'Al_A', 'type': 'Al', 'coords': np.array([0.3182, 0.2158, 0.2500])},
        {'label': 'Al_B', 'type': 'Al', 'coords': np.array([0.0000, 0.3662, 0.1030])},
        {'label': 'Al_C', 'type': 'Al', 'coords': np.array([0.1743, 0.0000, 0.0000])},
        {'label': 'Re_A', 'type': 'Re', 'coords': np.array([0.0000, 0.0445, 0.2500])},
    ]

    # Helper function to apply a symmetry operation
    def apply_sym_op(coords, op_str):
        x, y, z = coords[0], coords[1], coords[2]
        # Using eval is generally unsafe, but acceptable here as the operation strings are fixed.
        op_parts = op_str.split(',')
        new_x = eval(op_parts[0], {"x": x, "y": y, "z": z})
        new_y = eval(op_parts[1], {"x": x, "y": y, "z": z})
        new_z = eval(op_parts[2], {"x": x, "y": y, "z": z})
        return np.array([new_x, new_y, new_z])

    # Part 2: Generate all atoms in the unit cell
    all_atoms_in_cell = []
    seen_coords = set()
    for atom in au_atoms:
        for op in sym_ops_str:
            new_coords_raw = apply_sym_op(atom['coords'], op)
            # Normalize coordinates to be inside the unit cell [0, 1)
            new_coords_norm = new_coords_raw - np.floor(new_coords_raw)
            # Round to 4 decimal places to avoid float precision issues when checking for duplicates
            coord_tuple = tuple(np.round(new_coords_norm, 4))
            if coord_tuple not in seen_coords:
                seen_coords.add(coord_tuple)
                all_atoms_in_cell.append({'type': atom['type'], 'coords': new_coords_norm})

    # Part 3: Find coordination environment for each unique atom
    results = []
    for center_atom in au_atoms:
        center_coords = center_atom['coords']
        
        neighbors = []
        # Check for neighbors in a 3x3x3 supercell to handle cell boundaries
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    translation_vector = np.array([i, j, k])
                    for atom in all_atoms_in_cell:
                        # Calculate distance in Cartesian coordinates
                        delta_frac = center_coords - (atom['coords'] + translation_vector)
                        delta_cart = delta_frac * cell_lengths
                        distance = np.linalg.norm(delta_cart)

                        # Store if it's a neighbor (not the atom itself) within a reasonable cutoff
                        if 0.01 < distance < 4.0:
                             neighbors.append({'type': atom['type'], 'dist': distance})
        
        # Part 4: Identify the first coordination shell
        if not neighbors:
            continue
        
        neighbors.sort(key=lambda n: n['dist'])
        
        first_shell_neighbors = []
        if neighbors:
            first_shell_neighbors.append(neighbors[0])
            # A distance gap is defined as a jump of more than 15%
            for i in range(1, len(neighbors)):
                if neighbors[i]['dist'] / neighbors[i-1]['dist'] < 1.15:
                    first_shell_neighbors.append(neighbors[i])
                else:
                    break  # Found the first significant gap

        # Part 5: Formulate the result
        shell_counts = Counter(n['type'] for n in first_shell_neighbors)
        
        # Sort elements alphabetically (Al, Re) for a consistent formula
        formula_parts = []
        for element in sorted(shell_counts.keys()):
            count = shell_counts[element]
            formula_parts.append(f"{element}{count}" if count > 1 else element)
        formula = "".join(formula_parts)
        
        results.append(f"{center_atom['label']}, {formula}")
        
    print("; ".join(results))

if __name__ == "__main__":
    find_coordination_polyhedra()