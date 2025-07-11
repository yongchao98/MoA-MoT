import numpy as np
from collections import Counter

def solve():
    """
    This function determines the coordination polyhedra for the ReAl12 crystal structure.
    It follows these steps:
    1.  Define crystal structure data: lattice parameters, asymmetric atom positions, and symmetry operations.
    2.  Generate all unique atom positions in a single unit cell by applying symmetry operations to the asymmetric atoms.
    3.  Create a 3x3x3 supercell by translating the unit cell atoms to account for periodic boundary conditions. This ensures all nearest neighbors are found.
    4.  For each unique atom in the asymmetric unit (the "central atom"):
        a. Calculate the distances to all other atoms in the supercell.
        b. Sort these distances to identify nearest neighbors.
        c. The first coordination shell is determined by finding a significant gap in the sorted distances. A gap is identified when the ratio of two consecutive distances is greater than a threshold (1.2).
        d. The elements of the atoms in this first shell are counted.
        e. A chemical formula for the coordination polyhedron is constructed from these counts.
    5.  The results for each unique central atom type are collected.
    6.  The final, unique polyhedra are printed in the specified format: "central atom, formula; ...".
    """
    
    # 1. Define crystal structure data
    lattice = {
        'a': 7.6091, 'b': 6.6117, 'c': 9.0230,
        'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0
    }
    
    sym_ops_str = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2',
        '-x+1/2,y+1/2,z', '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2',
        'x+1/2,y+1/2,z', 'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]
    
    asym_atoms = [
        {'label': 'Al_A', 'symbol': 'Al', 'coords': np.array([0.3182, 0.2158, 0.2500])},
        {'label': 'Al_B', 'symbol': 'Al', 'coords': np.array([0.0000, 0.3662, 0.1030])},
        {'label': 'Al_C', 'symbol': 'Al', 'coords': np.array([0.1743, 0.0000, 0.0000])},
        {'label': 'Re_A', 'symbol': 'Re', 'coords': np.array([0.0000, 0.0445, 0.2500])},
    ]

    def get_sym_ops(ops_str_list):
        operations = []
        for op_str in ops_str_list:
            ops = op_str.replace(' ', '').split(',')
            # Use a closure to capture the correct op_x, op_y, op_z for each lambda
            op_func = lambda x, y, z, op_x=ops[0], op_y=ops[1], op_z=ops[2]: np.array([
                eval(op_x, {"x": x, "y": y, "z": z}),
                eval(op_y, {"x": x, "y": y, "z": z}),
                eval(op_z, {"x": x, "y": y, "z": z})
            ])
            operations.append(op_func)
        return operations

    sym_ops = get_sym_ops(sym_ops_str)

    # 2. Generate all unique atoms in the unit cell
    unit_cell_atoms = []
    seen_coords = []
    tolerance = 1e-4
    for atom in asym_atoms:
        for op in sym_ops:
            new_coords = op(atom['coords'][0], atom['coords'][1], atom['coords'][2])
            new_coords_mod = new_coords - np.floor(new_coords)
            
            is_new = True
            for seen_coord in seen_coords:
                delta = new_coords_mod - seen_coord
                delta -= np.round(delta) # periodic boundary condition
                if np.linalg.norm(delta) < tolerance:
                    is_new = False
                    break
            
            if is_new:
                unit_cell_atoms.append({'symbol': atom['symbol'], 'coords': new_coords_mod})
                seen_coords.append(new_coords_mod)

    # 3. Create a 3x3x3 supercell
    supercell_atoms = []
    for atom in unit_cell_atoms:
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    supercell_atoms.append({
                        'symbol': atom['symbol'],
                        'coords': atom['coords'] + np.array([i, j, k])
                    })

    # Transformation matrix from fractional to Cartesian coordinates (for orthorhombic)
    lat_matrix = np.diag([lattice['a'], lattice['b'], lattice['c']])
    
    results = []
    
    # 4. Iterate through each unique asymmetric atom
    for central_atom in asym_atoms:
        central_coords_frac = central_atom['coords']
        
        distances = []
        for neighbor_atom in supercell_atoms:
            delta_frac = neighbor_atom['coords'] - central_coords_frac
            delta_cart = np.dot(lat_matrix, delta_frac)
            dist = np.linalg.norm(delta_cart)

            if dist > tolerance:
                distances.append({'dist': dist, 'symbol': neighbor_atom['symbol']})

        distances.sort(key=lambda x: x['dist'])

        # 4c. Determine first coordination shell by finding a large gap
        first_shell_neighbors = []
        if distances:
            first_shell_neighbors.append(distances[0])
            for i in range(1, len(distances)):
                # If the ratio to the previous distance is large, it's a new shell
                if distances[i]['dist'] / distances[i-1]['dist'] > 1.20:
                    break
                first_shell_neighbors.append(distances[i])

        # 4d. Count elements in the shell
        neighbor_symbols = [n['symbol'] for n in first_shell_neighbors]
        counts = Counter(neighbor_symbols)
        
        # 4e. Create the formula string
        formula_parts = []
        # Define element order for consistent formula representation
        element_order = ['Re', 'Al']
        for elem in element_order:
            if elem in counts:
                count = counts[elem]
                if count > 1:
                    formula_parts.append(f"{elem}{count}")
                elif count == 1:
                    formula_parts.append(elem)

        formula = "".join(formula_parts)
        
        # 5. Add to results
        results.append(f"{central_atom['symbol']}, {formula}")

    # 6. Print unique, sorted results
    unique_results = sorted(list(set(results)))
    print("; ".join(unique_results))

solve()