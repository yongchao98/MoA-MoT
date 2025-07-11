import numpy as np
from collections import Counter

def find_coordination_polyhedra():
    """
    Analyzes a crystal structure to find the coordination polyhedra for each unique atom.
    The script is configured with data for ReAl12.
    """
    
    # 1. Parse Crystal Data from the CIF file for ReAl12
    a, b, c = 7.6091, 6.6117, 9.0230
    alpha, beta, gamma = 90.0, 90.0, 90.0

    symop_strings = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2', '-x+1/2,y+1/2,z',
        '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2', 'x+1/2,y+1/2,z',
        'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]

    # Basis atoms from the asymmetric unit
    basis_atoms = [
        {'label': 'Al_A', 'type': 'Al', 'pos': np.array([0.3182, 0.2158, 0.2500])},
        {'label': 'Al_B', 'type': 'Al', 'pos': np.array([0.0000, 0.3662, 0.1030])},
        {'label': 'Al_C', 'type': 'Al', 'pos': np.array([0.1743, 0.0000, 0.0000])},
        {'label': 'Re_A', 'type': 'Re', 'pos': np.array([0.0000, 0.0445, 0.2500])}
    ]

    # Helper function to parse symmetry operation strings
    def parse_symop(s):
        x_op, y_op, z_op = s.split(',')
        return lambda x, y, z: (eval(x_op), eval(y_op), eval(z_op))

    symops = [parse_symop(s) for s in symop_strings]

    # Transformation matrix for fractional to Cartesian coordinates (orthorhombic case)
    trans_matrix = np.diag([a, b, c])

    # 2. Generate all atoms in the full unit cell
    def generate_full_cell():
        full_cell_atoms = []
        unique_pos_set = set()
        for atom_info in basis_atoms:
            pos = atom_info['pos']
            for op in symops:
                new_pos = op(pos[0], pos[1], pos[2])
                new_pos_norm = np.array([p - np.floor(p) for p in new_pos])
                
                # Check for duplicates using rounded positions
                rounded_pos_tuple = tuple(np.round(new_pos_norm, 4))
                if rounded_pos_tuple not in unique_pos_set:
                    unique_pos_set.add(rounded_pos_tuple)
                    full_cell_atoms.append({'type': atom_info['type'], 'pos': new_pos_norm})
        return full_cell_atoms

    all_atoms_in_cell = generate_full_cell()

    results = []
    # 3. Iterate through each unique basis atom as the polyhedron center
    for central_atom_info in basis_atoms:
        center_type = central_atom_info['type']
        center_pos = central_atom_info['pos']
        
        # 4. Calculate distances to all neighbors in a 3x3x3 supercell
        neighbors = []
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    offset = np.array([i, j, k])
                    for atom in all_atoms_in_cell:
                        diff_fract = atom['pos'] + offset - center_pos
                        dist = np.linalg.norm(trans_matrix @ diff_fract)
                        
                        if dist > 1e-4: # Exclude the atom itself
                            neighbors.append({'type': atom['type'], 'dist': dist})
        
        neighbors.sort(key=lambda n: n['dist'])
        
        # 5. Identify the first coordination shell
        if not neighbors:
            continue
            
        first_shell = []
        if neighbors:
            first_dist = neighbors[0]['dist']
            # Use a threshold factor to find atoms in the first shell
            cutoff_dist = first_dist * 1.25 
            for neighbor in neighbors:
                if neighbor['dist'] <= cutoff_dist:
                    first_shell.append(neighbor)
                else:
                    break
        
        # 6. Determine the polyhedron formula
        counts = Counter(n['type'] for n in first_shell)
        
        # Define element order for consistent formula string (e.g., ReAl not AlRe)
        element_order = ['Re', 'Al']
        sorted_elements = sorted(counts.keys(), key=lambda x: element_order.index(x) if x in element_order else len(element_order))
        
        formula_parts = []
        for elem in sorted_elements:
            num = counts[elem]
            formula_parts.append(f"{elem}{num if num > 1 else ''}")
        
        formula = "".join(formula_parts)
        results.append(f"{center_type}, {formula}")

    # Format the final output string, removing duplicates and sorting
    unique_results = sorted(list(set(results)))
    print("; ".join(unique_results))

find_coordination_polyhedra()