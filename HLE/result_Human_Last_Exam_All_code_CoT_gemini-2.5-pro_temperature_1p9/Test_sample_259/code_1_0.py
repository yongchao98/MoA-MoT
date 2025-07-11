import numpy as np
from collections import Counter

def find_coordination_polyhedra():
    """
    Determines the coordination polyhedra for the ReAl6 crystal structure based on provided CIF data.
    The script parses crystal information, generates all atomic positions in the unit cell,
    calculates interatomic distances considering periodic boundary conditions, identifies
    the first coordination shell using a gap-detection method, and formats the resulting
    polyhedra as requested.
    """
    # Step 1: Define Crystal Information from the CIF data
    a, b, c = 7.609100, 6.611700, 9.023000
    # For an orthorhombic cell, the lattice matrix is diagonal
    lat_matrix = np.array([[a, 0, 0], [0, b, 0], [0, 0, c]])

    # Define the atoms in the asymmetric unit
    asymmetric_atoms = [
        {'label': 'Al_A', 'type': 'Al', 'coords': np.array([0.3182, 0.2158, 0.2500])},
        {'label': 'Al_B', 'type': 'Al', 'coords': np.array([0.0000, 0.3662, 0.1030])},
        {'label': 'Al_C', 'type': 'Al', 'coords': np.array([0.1743, 0.0000, 0.0000])},
        {'label': 'Re_A', 'type': 'Re', 'coords': np.array([0.0000, 0.0445, 0.2500])},
    ]

    # Define the symmetry operations
    sym_ops_str = [
        "-x+1/2,-y+1/2,-z", "-x+1/2,-y+1/2,z+1/2", "-x+1/2,y+1/2,-z+1/2",
        "-x+1/2,y+1/2,z", "-x,-y,-z", "-x,-y,z+1/2",
        "-x,y,-z+1/2", "-x,y,z", "x+1/2,-y+1/2,-z",
        "x+1/2,-y+1/2,z+1/2", "x+1/2,y+1/2,-z+1/2", "x+1/2,y+1/2,z",
        "x,-y,-z", "x,-y,z+1/2", "x,y,-z+1/2", "x,y,z"
    ]

    def apply_sym_op(coords, op_str):
        x, y, z = coords
        # Use a controlled scope for eval to safely parse symmetry operation strings
        vars_map = {'x': x, 'y': y, 'z': z}
        op_parts = op_str.split(',')
        new_coords = [eval(part, {"__builtins__": None}, vars_map) for part in op_parts]
        return np.array(new_coords)

    # Step 2: Generate all atom positions in the unit cell
    all_atoms = []
    seen_coords_tuples = set()
    for atom in asymmetric_atoms:
        for op_str in sym_ops_str:
            new_coords_raw = apply_sym_op(atom['coords'], op_str)
            # Use modulo and rounding to identify unique atomic positions
            new_coords_mod = new_coords_raw % 1.0
            new_coords_tuple = tuple(np.round(new_coords_mod, decimals=4))
            
            if new_coords_tuple not in seen_coords_tuples:
                seen_coords_tuples.add(new_coords_tuple)
                all_atoms.append({'type': atom['type'], 'coords': new_coords_raw})

    # Steps 3, 4, 5: Find coordination for each unique central atom
    results = []
    # Sort to ensure consistent output order
    central_atoms = sorted(asymmetric_atoms, key=lambda x: x['label'])

    for central_atom in central_atoms:
        c_coords_frac = central_atom['coords']
        
        # Calculate distances to all other atoms in a 3x3x3 supercell
        neighbors = []
        for atom in all_atoms:
            for i in range(-1, 2):
                for j in range(-1, 2):
                    for k in range(-1, 2):
                        n_coords_frac = atom['coords'] + np.array([i, j, k])
                        delta_frac = n_coords_frac - c_coords_frac
                        delta_cart = np.dot(delta_frac, lat_matrix)
                        dist = np.linalg.norm(delta_cart)

                        # Exclude distance to self
                        if dist > 0.01:
                            neighbors.append({'dist': dist, 'type': atom['type']})
        
        neighbors.sort(key=lambda x: x['dist'])

        if not neighbors:
            continue
            
        unique_sorted_distances = sorted(list(set(round(n['dist'], 4) for n in neighbors)))
        
        # Determine cutoff by finding the largest gap in distances (up to 4.0 Å)
        cutoff = unique_sorted_distances[0] * 1.4  # Fallback cutoff
        relevant_distances = [d for d in unique_sorted_distances if d < 4.0]

        if len(relevant_distances) > 1:
            gaps = [relevant_distances[i] - relevant_distances[i-1] for i in range(1, len(relevant_distances))]
            if gaps:
                max_gap = max(gaps)
                if max_gap > 0.4:  # A gap > 0.4 Å is considered significant
                    max_gap_index = gaps.index(max_gap)
                    # Set cutoff in the middle of the largest gap
                    cutoff = relevant_distances[max_gap_index] + max_gap / 2.0
        
        first_shell_neighbors = [n for n in neighbors if n['dist'] <= cutoff]
        
        # Count atom types in the coordination shell to create the formula
        counts = Counter(n['type'] for n in first_shell_neighbors)
        # Sort alphabetically by element symbol for a standard formula
        formula_parts = sorted([f"{el}{count}" for el, count in counts.items()])
        formula = "".join(formula_parts)
        
        results.append(f"{central_atom['type']}, {formula}")

    # Step 6: Format and print the final result
    print("; ".join(results))

find_coordination_polyhedra()