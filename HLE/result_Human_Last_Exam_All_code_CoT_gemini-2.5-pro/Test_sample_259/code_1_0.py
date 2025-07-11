import numpy as np
from collections import Counter

def solve_coordination():
    """
    Analyzes a crystal structure from CIF data to find coordination polyhedra.
    """
    # 1. Parse Crystal Data from the problem description
    a, b, c = 7.6091, 6.6117, 9.0230
    lattice_vectors = np.array([a, b, c])

    # Unique atom sites in the asymmetric unit
    base_atoms = [
        {'label': 'Al_A', 'symbol': 'Al', 'coords': np.array([0.3182, 0.2158, 0.2500])},
        {'label': 'Al_B', 'symbol': 'Al', 'coords': np.array([0.0000, 0.3662, 0.1030])},
        {'label': 'Al_C', 'symbol': 'Al', 'coords': np.array([0.1743, 0.0000, 0.0000])},
        {'label': 'Re_A', 'symbol': 'Re', 'coords': np.array([0.0000, 0.0445, 0.2500])},
    ]

    # Symmetry operations for space group Cmcm (No. 63)
    sym_ops_str = [
        "-x+1/2,-y+1/2,-z", "-x+1/2,-y+1/2,z+1/2", "-x+1/2,y+1/2,-z+1/2", "-x+1/2,y+1/2,z",
        "-x,-y,-z", "-x,-y,z+1/2", "-x,y,-z+1/2", "-x,y,z",
        "x+1/2,-y+1/2,-z", "x+1/2,-y+1/2,z+1/2", "x+1/2,y+1/2,-z+1/2", "x+1/2,y+1/2,z",
        "x,-y,-z", "x,-y,z+1/2", "x,y,-z+1/2", "x,y,z",
    ]

    def parse_sym_op(op_str):
        """Parses a symmetry operation string into a callable function."""
        parts = op_str.split(',')
        op_lambda = lambda x, y, z: (
            eval(parts[0], {"x": x, "y": y, "z": z}),
            eval(parts[1], {"x": x, "y": y, "z": z}),
            eval(parts[2], {"x": x, "y": y, "z": z})
        )
        return op_lambda

    sym_ops = [parse_sym_op(s) for s in sym_ops_str]

    # 2. Generate Full Unit Cell
    unit_cell_atoms = []
    seen_coords_tuples = set()
    TOL = 1e-4
    for atom in base_atoms:
        for op in sym_ops:
            x, y, z = atom['coords']
            new_coords_raw = op(x, y, z)
            new_coords = np.array([coord % 1.0 for coord in new_coords_raw])
            
            # Use a rounded tuple for checking existence in the set
            new_coords_tuple = tuple(np.round(new_coords, 4))
            
            # Check for duplicates using minimum image distance to be robust
            is_duplicate = False
            for seen_tuple in seen_coords_tuples:
                diff = new_coords - np.array(seen_tuple)
                diff -= np.round(diff) # minimum image convention
                if np.linalg.norm(diff) < TOL:
                    is_duplicate = True
                    break
            
            if not is_duplicate:
                unit_cell_atoms.append({'symbol': atom['symbol'], 'coords': new_coords})
                seen_coords_tuples.add(new_coords_tuple)

    # 3. Create Supercell
    supercell_atoms = []
    for atom in unit_cell_atoms:
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    translation = np.array([i, j, k])
                    supercell_atoms.append({
                        'symbol': atom['symbol'],
                        'coords': atom['coords'] + translation
                    })

    # 4-7. Find and describe polyhedra for each unique atom
    final_results = []
    for central_atom in base_atoms:
        distances = []
        central_coords = central_atom['coords']
        
        for neighbor in supercell_atoms:
            diff_frac = central_coords - neighbor['coords']
            dist = np.linalg.norm(diff_frac * lattice_vectors)
            
            if dist > TOL:  # Exclude the atom itself
                distances.append({'dist': dist, 'symbol': neighbor['symbol']})

        distances.sort(key=lambda x: x['dist'])
        
        # Identify coordination shell by finding a large gap in distances
        if not distances: continue
        shell_cutoff = distances[0]['dist'] * 1.25 # Heuristic cutoff
        sorted_dists = [d['dist'] for d in distances]
        for i in range(len(sorted_dists) - 1):
            if sorted_dists[i+1] > sorted_dists[i] * 1.25:
                shell_cutoff = sorted_dists[i]
                break
        
        first_shell_symbols = [d['symbol'] for d in distances if d['dist'] <= shell_cutoff + TOL]
        
        # Determine polyhedron formula
        counts = Counter(first_shell_symbols)
        formula_parts = []
        for symbol in sorted(counts.keys()):
            count = counts[symbol]
            formula_parts.append(f"{symbol}{count if count > 1 else ''}")
        
        formula = "".join(formula_parts)
        final_results.append(f"{central_atom['symbol']}, {formula}")
    
    # 8. Format and print the final output
    print("; ".join(final_results))

solve_coordination()