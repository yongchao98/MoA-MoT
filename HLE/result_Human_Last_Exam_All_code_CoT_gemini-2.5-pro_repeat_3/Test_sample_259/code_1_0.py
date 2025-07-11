import numpy as np
from collections import Counter

def solve_coordination_polyhedra():
    """
    Analyzes a crystal structure to determine the coordination polyhedra for each unique atom.
    """
    # 1. Define crystal structure data from the CIF file
    cell_params = [7.6091, 6.6117, 9.0230]  # a, b, c for orthorhombic cell
    unique_atoms = [
        ("Al_A", "Al", np.array([0.3182, 0.2158, 0.2500])),
        ("Al_B", "Al", np.array([0.0000, 0.3662, 0.1030])),
        ("Al_C", "Al", np.array([0.1743, 0.0000, 0.0000])),
        ("Re_A", "Re", np.array([0.0000, 0.0445, 0.2500]))
    ]
    sym_ops_str = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2',
        '-x+1/2,y+1/2,z', '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2',
        '-x,y,z', 'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2',
        'x+1/2,y+1/2,-z+1/2', 'x+1/2,y+1/2,z', 'x,-y,-z',
        'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]

    def apply_sym_op(op_str, pos):
        x, y, z = pos[0], pos[1], pos[2]
        new_x = eval(op_str.split(',')[0], {"x": x, "y": y, "z": z})
        new_y = eval(op_str.split(',')[1], {"x": x, "y": y, "z": z})
        new_z = eval(op_str.split(',')[2], {"x": x, "y": y, "z": z})
        return np.array([new_x, new_y, new_z])

    def get_dist(p1, p2, cell):
        delta_frac = p1 - p2
        delta_cart = delta_frac * np.array(cell)
        return np.linalg.norm(delta_cart)

    # 2. Generate all atoms in the base unit cell
    base_cell_atoms = []
    seen_coords = set()
    for _, symbol, pos in unique_atoms:
        for op_str in sym_ops_str:
            new_pos = apply_sym_op(op_str, pos)
            norm_pos = new_pos - np.floor(new_pos)
            coord_tuple = tuple(np.round(norm_pos, 4))
            if coord_tuple not in seen_coords:
                seen_coords.add(coord_tuple)
                base_cell_atoms.append({'symbol': symbol, 'pos': new_pos})

    # 3. Generate a 3x3x3 supercell
    supercell_atoms = []
    for atom in base_cell_atoms:
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    supercell_atoms.append({
                        'symbol': atom['symbol'],
                        'pos': atom['pos'] + np.array([i, j, k])
                    })

    # 4. Analyze coordination for each unique central atom
    polyhedra = {}
    for label, center_symbol, center_pos in unique_atoms:
        # 5. Calculate distances to all neighbors
        neighbors = []
        for neighbor in supercell_atoms:
            dist = get_dist(center_pos, neighbor['pos'], cell_params)
            if dist > 0.01:  # Exclude the atom itself
                neighbors.append({'dist': dist, 'symbol': neighbor['symbol']})
        
        neighbors.sort(key=lambda n: n['dist'])

        if not neighbors:
            continue

        # 6. Find the first coordination shell by finding a significant distance gap
        cutoff_dist = neighbors[0]['dist'] * 1.25 # Initial guess
        for i in range(len(neighbors) - 1):
            if neighbors[i+1]['dist'] > neighbors[i]['dist'] * 1.2:
                cutoff_dist = (neighbors[i]['dist'] + neighbors[i+1]['dist']) / 2
                break
        
        first_shell_neighbors = [n for n in neighbors if n['dist'] <= cutoff_dist]
        
        # 7. Count elements and generate the formula
        shell_counts = Counter(n['symbol'] for n in first_shell_neighbors)
        formula_parts = []
        for elem, count in sorted(shell_counts.items()):
            formula_parts.append(f"{elem}{count if count > 1 else ''}")
        formula = "".join(formula_parts)
        
        # Store unique polyhedra descriptions, keyed by the original site label
        polyhedra[label] = f"{center_symbol}, {formula}"
    
    # 8. Format the final output string, ensuring unique environments are listed
    # Use a dictionary to get unique values, then sort for consistent output
    final_results = sorted(list(set(polyhedra.values())))
    print("; ".join(final_results))

solve_coordination_polyhedra()