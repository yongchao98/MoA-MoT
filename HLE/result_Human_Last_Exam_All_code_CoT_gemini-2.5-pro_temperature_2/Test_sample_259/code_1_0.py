import math
from collections import Counter

def find_coordination_polyhedra():
    """
    This script determines the coordination polyhedra for a given crystal structure.
    
    1.  Define crystal structure data (cell, atoms, symmetry) from the CIF input.
    2.  Implement a function to apply symmetry operations to generate the full unit cell.
    3.  Create a 3x3x3 supercell by translating the unit cell atoms.
    4.  For each crystallographically unique atom (Re_A, Al_A, etc.):
        a. Calculate the distances to all other atoms in the supercell.
        b. Sort the neighbors by distance.
        c. Determine the first coordination shell by finding the largest gap in distances.
        d. Count atom types within the first shell to get the polyhedron's formula.
    5.  Format the results as "central_atom, formula;" for each site and print the final string.
    """
    # 1. Define crystal structure data
    cell_lengths = {'a': 7.6091, 'b': 6.6117, 'c': 9.0230}
    
    sym_ops_str = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2', '-x+1/2,y+1/2,z',
        '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2', 'x+1/2,y+1/2,z',
        'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]

    unique_atoms = [
        {'label': 'Al_A', 'symbol': 'Al', 'coords': [0.3182, 0.2158, 0.2500]},
        {'label': 'Al_B', 'symbol': 'Al', 'coords': [0.0000, 0.3662, 0.1030]},
        {'label': 'Al_C', 'symbol': 'Al', 'coords': [0.1743, 0.0000, 0.0000]},
        {'label': 'Re_A', 'symbol': 'Re', 'coords': [0.0000, 0.0445, 0.2500]}
    ]

    # 2. Generate full unit cell
    def apply_sym_op(coords, op_str):
        x, y, z = coords
        op_x, op_y, op_z = op_str.split(',')
        # Safely evaluate the expression
        new_x = eval(op_x, {'x': x, 'y': y, 'z': z})
        new_y = eval(op_y, {'x': x, 'y': y, 'z': z})
        new_z = eval(op_z, {'x': x, 'y': y, 'z': z})
        return [new_x, new_y, new_z]

    unit_cell_atoms = []
    seen_coords = set()
    for atom_info in unique_atoms:
        for op in sym_ops_str:
            new_coords = apply_sym_op(atom_info['coords'], op)
            norm_coords = tuple(round(c - math.floor(c) + 1e-9, 5) % 1.0 for c in new_coords)
            if norm_coords not in seen_coords:
                seen_coords.add(norm_coords)
                unit_cell_atoms.append({'symbol': atom_info['symbol'], 'coords': list(norm_coords)})
    
    # 3. Construct a supercell
    supercell_atoms = []
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                for atom in unit_cell_atoms:
                    supercell_atoms.append({
                        'symbol': atom['symbol'],
                        'coords': [atom['coords'][0] + i, atom['coords'][1] + j, atom['coords'][2] + k]
                    })
    
    def get_distance(p1, p2, cell):
        dx = (p1[0] - p2[0]) * cell['a']
        dy = (p1[1] - p2[1]) * cell['b']
        dz = (p1[2] - p2[2]) * cell['c']
        return math.sqrt(dx**2 + dy**2 + dz**2)

    results = []
    # 4. Find coordination for each unique atom
    for central_atom_info in unique_atoms:
        center_coords = central_atom_info['coords']
        
        # a. Calculate distances
        distances = []
        for neighbor in supercell_atoms:
            dist = get_distance(center_coords, neighbor['coords'], cell_lengths)
            if dist > 1e-6:
                distances.append({'symbol': neighbor['symbol'], 'dist': dist})
        
        # b. Sort by distance
        distances.sort(key=lambda x: x['dist'])

        # c. Find coordination shell by max gap method
        first_shell_atoms = []
        if distances:
            gaps = [distances[i+1]['dist'] - distances[i]['dist'] for i in range(len(distances)-1)]
            if not gaps:
                 cn = len(distances)
            else:
                max_gap = 0
                max_gap_index = -1
                for i, gap in enumerate(gaps):
                    # Find a significant gap
                    if gap > 0.5: # A gap of 0.5 Angstrom is significant
                        max_gap_index = i
                        break
                if max_gap_index != -1:
                    cn = max_gap_index + 1
                else: # No large gap, may be a complex case, use ratio
                    cn = 0
                    for i in range(len(distances) - 1):
                        cn = i + 1
                        if distances[i+1]['dist'] > distances[i]['dist'] * 1.25:
                           break
                    else:
                        cn += 1

            first_shell_atoms = distances[:cn]

        # d. & e. Get formula
        neighbor_counts = Counter(atom['symbol'] for atom in first_shell_atoms)
        formula_parts = []
        for symbol in sorted(neighbor_counts.keys()):
            count = neighbor_counts[symbol]
            formula_parts.append(f"{symbol}{count if count > 1 else ''}")
        
        formula = "".join(formula_parts)
        central_atom_symbol = central_atom_info['label'].split('_')[0]
        results.append(f"{central_atom_symbol}, {formula}")

    # 6. Format and print output
    # Since Al_B and Al_C might be equivalent by relabeling axes but are crystallographically
    # distinct as given, we present results for all four unique sites.
    # Joining with '; ' and adding a trailing one if requested.
    final_output = "; ".join(results) + ";"
    print(final_output)

find_coordination_polyhedra()