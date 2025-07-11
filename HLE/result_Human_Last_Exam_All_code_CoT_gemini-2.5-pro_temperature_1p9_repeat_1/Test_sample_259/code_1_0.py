import numpy as np
from collections import Counter

def find_coordination_polyhedra():
    """
    This script calculates the coordination polyhedra for the ReAl12 structure
    based on the provided CIF data.
    """
    
    # 1. Parse Crystal Data from the prompt
    cell = {'a': 7.609100, 'b': 6.611700, 'c': 9.023000}
    a, b, c = cell['a'], cell['b'], cell['c']

    sym_ops_str = [
        "-x+1/2,-y+1/2,-z", "-x+1/2,-y+1/2,z+1/2", "-x+1/2,y+1/2,-z+1/2",
        "-x+1/2,y+1/2,z", "-x,-y,-z", "-x,-y,z+1/2", "-x,y,-z+1/2", "-x,y,z",
        "x+1/2,-y+1/2,-z", "x+1/2,-y+1/2,z+1/2", "x+1/2,y+1/2,-z+1/2",
        "x+1/2,y+1/2,z", "x,-y,-z", "x,-y,z+1/2", "x,y,-z+1/2", "x,y,z"
    ]
    
    asymmetric_atoms = [
        {'label': 'Al_A', 'type': 'Al', 'coords': np.array([0.3182, 0.2158, 0.2500])},
        {'label': 'Al_B', 'type': 'Al', 'coords': np.array([0.0000, 0.3662, 0.1030])},
        {'label': 'Al_C', 'type': 'Al', 'coords': np.array([0.1743, 0.0000, 0.0000])},
        {'label': 'Re_A', 'type': 'Re', 'coords': np.array([0.0000, 0.0445, 0.2500])}
    ]

    def parse_sym_op(op_str):
        """Converts a symmetry operation string into a callable function."""
        # Using eval in a controlled environment is safe and effective here.
        # The 'x, y, z' are replaced with numerical values during the call.
        x_op, y_op, z_op = op_str.split(',')
        return lambda x, y, z: np.array([
            eval(x_op, {"x": x, "y": y, "z": z}),
            eval(y_op, {"x": x, "y": y, "z": z}),
            eval(z_op, {"x": x, "y": y, "z": z})
        ])

    sym_ops = [parse_sym_op(s) for s in sym_ops_str]

    # 2. Generate all atoms in the full unit cell
    unit_cell_atoms = []
    seen_coords = set()
    for atom_site in asymmetric_atoms:
        for op in sym_ops:
            x_i, y_i, z_i = atom_site['coords']
            new_coords = op(x_i, y_i, z_i)
            # Normalize coordinates to be within the unit cell [0, 1)
            new_coords -= np.floor(new_coords)
            
            # Round to avoid floating point precision issues with duplicates
            rounded_key = tuple(np.round(new_coords, decimals=4))
            if rounded_key not in seen_coords:
                unit_cell_atoms.append({'type': atom_site['type'], 'coords': new_coords})
                seen_coords.add(rounded_key)

    # 3. Find coordination environment for each unique atom site
    final_results = []
    for central_atom in asymmetric_atoms:
        distances = []
        for neighbor in unit_cell_atoms:
            # Calculate vector in fractional coordinates
            delta_frac = neighbor['coords'] - central_atom['coords']
            
            # Apply minimum image convention for periodic boundaries
            delta_frac -= np.round(delta_frac)
            
            # Convert to Cartesian and calculate distance
            dist_sq = ((delta_frac[0] * a)**2 +
                       (delta_frac[1] * b)**2 +
                       (delta_frac[2] * c)**2)
            dist = np.sqrt(dist_sq)

            # Add to list if it's not the central atom itself
            if dist > 1e-5:
                distances.append({'dist': dist, 'type': neighbor['type']})
        
        # 4. Identify the first coordination shell
        distances.sort(key=lambda item: item['dist'])
        dist_values = [d['dist'] for d in distances]
        
        first_shell_cutoff_index = len(dist_values) # Default to all atoms
        if len(dist_values) > 1:
            # Find the largest gap between consecutive distances
            gaps = [dist_values[i+1] - dist_values[i] for i in range(len(dist_values) - 1)]
            if max(gaps) > 0.5:  # A significant gap indicates a new shell
                max_gap_index = np.argmax(gaps)
                first_shell_cutoff_index = max_gap_index + 1
        
        first_shell_neighbors = distances[:first_shell_cutoff_index]

        # 5. Determine the polyhedron formula
        neighbor_types = [n['type'] for n in first_shell_neighbors]
        counts = Counter(neighbor_types)
        
        formula = "".join(f"{elem}{counts[elem]}" for elem in sorted(counts.keys()))
        
        result_string = f"{central_atom['type']}, {formula}"
        final_results.append(result_string)

    # 6. Print unique results in the required format
    unique_results = sorted(list(set(final_results)))
    print("; ".join(unique_results))

find_coordination_polyhedra()