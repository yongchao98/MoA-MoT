import numpy as np
from collections import Counter

def solve_coordination():
    """
    Analyzes a crystal structure from CIF data to determine coordination polyhedra.
    """
    # 1. Parse CIF Data
    cell_lengths = {'a': 7.609100, 'b': 6.611700, 'c': 9.023000}
    # Orthorhombic cell, angles are 90 degrees
    cell_angles_rad = {'alpha': np.pi/2, 'beta': np.pi/2, 'gamma': np.pi/2}

    asym_atoms_data = [
        {'label': 'Al_A', 'symbol': 'Al', 'coords': np.array([0.3182, 0.2158, 0.2500])},
        {'label': 'Al_B', 'symbol': 'Al', 'coords': np.array([0.0000, 0.3662, 0.1030])},
        {'label': 'Al_C', 'symbol': 'Al', 'coords': np.array([0.1743, 0.0000, 0.0000])},
        {'label': 'Re_A', 'symbol': 'Re', 'coords': np.array([0.0000, 0.0445, 0.2500])},
    ]

    sym_ops_str = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2', '-x+1/2,y+1/2,z',
        '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2', 'x+1/2,y+1/2,z',
        'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z',
    ]

    # Function to convert string operations to executable functions
    def parse_sym_op(op_str):
        return eval(f"lambda x, y, z: ({op_str})")

    sym_ops = [parse_sym_op(s) for s in sym_ops_str]

    # 2. Generate Full Unit Cell
    all_atoms_in_cell = []
    unique_positions = set()
    for atom in asym_atoms_data:
        for op in sym_ops:
            fx, fy, fz = op(*atom['coords'])
            # Normalize coordinates to be within [0, 1) to handle periodicity
            pos_norm = np.array([fx % 1, fy % 1, fz % 1])
            # Round for uniqueness check to handle floating point inaccuracies
            rounded_pos = tuple(np.round(pos_norm, 4))
            unique_key = (atom['symbol'], rounded_pos)
            
            if unique_key not in unique_positions:
                unique_positions.add(unique_key)
                all_atoms_in_cell.append({'symbol': atom['symbol'], 'frac_coords': pos_norm})

    final_results = []
    
    # Process each unique atom type from the asymmetric unit
    for central_atom in asym_atoms_data:
        # 3. Calculate Distances
        distances = []
        central_frac = central_atom['coords']
        
        # Iterate over a 3x3x3 supercell to handle periodic boundaries
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    offset = np.array([i, j, k])
                    for neighbor in all_atoms_in_cell:
                        delta_frac = neighbor['frac_coords'] + offset - central_frac
                        # For orthorhombic cell, Cartesian conversion is direct
                        delta_cart = delta_frac * np.array([cell_lengths['a'], cell_lengths['b'], cell_lengths['c']])
                        dist = np.linalg.norm(delta_cart)
                        
                        # Exclude the distance to the atom itself
                        if dist > 1e-4:
                            distances.append({'dist': dist, 'symbol': neighbor['symbol']})

        # Sort neighbors by distance
        distances.sort(key=lambda x: x['dist'])

        if not distances:
            continue

        # 4. Identify Coordination Shell by finding the largest gap
        max_ratio = 0
        cutoff_index = 0
        # Search for a gap within the first 25 neighbors
        search_range = min(len(distances) - 1, 25)
        for i in range(search_range):
            ratio = distances[i+1]['dist'] / distances[i]['dist']
            if ratio > max_ratio:
                max_ratio = ratio
                cutoff_index = i

        first_shell = distances[:cutoff_index + 1]

        # 5. Determine Polyhedron Formula
        neighbor_symbols = [d['symbol'] for d in first_shell]
        counts = Counter(neighbor_symbols)
        
        formula_parts = []
        # Sort symbols alphabetically for a consistent formula (Al, then Re)
        for symbol in sorted(counts.keys()):
            num = counts[symbol]
            formula_parts.append(f"{symbol}{num if num > 1 else ''}")
        formula = "".join(formula_parts)
        
        final_results.append(f"{central_atom['label']}, {formula}")
        
    # 6. Format and Print Output
    print("; ".join(final_results))

solve_coordination()