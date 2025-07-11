import numpy as np
from collections import Counter

def find_coordination_polyhedra():
    """
    This script determines the coordination polyhedra for each unique atom
    in the provided ReAl6 crystal structure.
    It calculates the nearest neighbors for each crystallographically unique
    atom and reports the resulting chemical formula of the coordination shell.
    """
    # Step 1: Define the crystal structure data from the CIF file.
    # Cell parameters for the orthorhombic cell
    a, b, c = 7.6091, 6.6117, 9.0230
    cell_vectors = np.array([a, b, c])

    # Asymmetric unit atoms (unique crystallographic sites)
    # Each atom is defined by its element symbol and fractional coordinates.
    asym_atoms = [
        {'symbol': 'Al', 'pos': np.array([0.3182, 0.2158, 0.2500])}, # Al_A
        {'symbol': 'Al', 'pos': np.array([0.0000, 0.3662, 0.1030])}, # Al_B
        {'symbol': 'Al', 'pos': np.array([0.1743, 0.0000, 0.0000])}, # Al_C
        {'symbol': 'Re', 'pos': np.array([0.0000, 0.0445, 0.2500])}  # Re_A
    ]

    # Symmetry operations for space group Cmcm (No. 63)
    sym_ops_str = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2', '-x+1/2,y+1/2,z',
        '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2', 'x+1/2,y+1/2,z',
        'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]

    # Helper function to apply a symmetry operation
    def apply_sym_op(op_str, coord):
        x, y, z = coord[0], coord[1], coord[2]
        # Use a safe evaluation method
        local_scope = {'x': x, 'y': y, 'z': z}
        parts = op_str.split(',')
        new_x = eval(parts[0], {"__builtins__": None}, local_scope)
        new_y = eval(parts[1], {"__builtins__": None}, local_scope)
        new_z = eval(parts[2], {"__builtins__": None}, local_scope)
        return np.array([new_x, new_y, new_z])

    # Step 2: Generate all atoms in the unit cell by applying symmetry operations.
    full_cell_atoms = []
    seen_coords = set()
    for atom in asym_atoms:
        for op_str in sym_ops_str:
            new_pos = apply_sym_op(op_str, atom['pos'])
            # Normalize coordinates to [0, 1) to identify unique positions
            norm_pos = new_pos % 1.0
            coord_tuple = tuple(np.round(norm_pos, 4))
            if coord_tuple not in seen_coords:
                seen_coords.add(coord_tuple)
                # Store the non-normalized position for distance calculations
                full_cell_atoms.append({'symbol': atom['symbol'], 'pos': new_pos})

    # Step 3: Create a 3x3x3 supercell to find all nearest neighbors.
    supercell_atoms = []
    for atom in full_cell_atoms:
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    supercell_atoms.append({
                        'symbol': atom['symbol'],
                        'pos': atom['pos'] + np.array([i, j, k])
                    })

    # Step 4: For each unique atom, find its coordination environment.
    final_results = []
    GAP_THRESHOLD = 0.5  # Angstroms, used to find the end of a coordination shell

    for central_atom in asym_atoms:
        central_pos = central_atom['pos']
        central_symbol = central_atom['symbol']

        # Calculate distances to all atoms in the supercell
        distances = []
        for neighbor_atom in supercell_atoms:
            # Skip the distance to the atom itself
            if np.allclose(central_pos, neighbor_atom['pos']):
                continue
            
            delta_frac = neighbor_atom['pos'] - central_pos
            delta_cart = delta_frac * cell_vectors
            dist = np.linalg.norm(delta_cart)
            
            if dist > 1e-4: # Avoid zero distance
                distances.append({'dist': dist, 'symbol': neighbor_atom['symbol']})

        # Sort neighbors by distance
        distances.sort(key=lambda x: x['dist'])

        # Identify the first coordination shell by finding a large gap in distances
        if not distances:
            continue
        
        first_shell = []
        last_dist = distances[0]['dist']
        first_shell.append(distances[0])

        for i in range(1, len(distances)):
            current_dist = distances[i]['dist']
            # If the distance to the previous atom is small, it's in the same shell
            if (current_dist - last_dist) < GAP_THRESHOLD:
                first_shell.append(distances[i])
            # Also include atoms at the exact same distance
            elif np.isclose(current_dist, last_dist):
                first_shell.append(distances[i])
            else:
                break  # A significant gap was found, the shell is complete
            last_dist = current_dist
            
        # Step 5: Count atoms in the shell and format the polyhedron formula.
        shell_symbols = [d['symbol'] for d in first_shell]
        counts = Counter(shell_symbols)
        
        formula_parts = []
        # Sort elements alphabetically for a consistent formula (e.g., Al before Re)
        for elem in sorted(counts.keys()):
            count = counts[elem]
            # Omit the number if the count is 1
            formula_parts.append(f"{elem}{count if count > 1 else ''}")
        
        formula = "".join(formula_parts)
        
        final_results.append(f"{central_symbol}, {formula}")

    # Step 6: Print the unique coordination environments found.
    # Using a set removes duplicate environments (e.g., if two Al sites have the same coordination)
    unique_results = sorted(list(set(final_results)))
    print("; ".join(unique_results))

find_coordination_polyhedra()