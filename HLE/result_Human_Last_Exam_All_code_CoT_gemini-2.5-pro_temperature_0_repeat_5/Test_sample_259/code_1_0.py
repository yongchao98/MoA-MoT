import numpy as np
from collections import Counter
import math

def solve():
    """
    This script determines the coordination polyhedra for the crystal structure of ReAl6.
    It follows the plan outlined above to calculate and report the coordination environment
    for each unique atom in the structure.
    """

    # 1. Crystal structure data from CIF
    cif_data = {
        "cell_lengths": np.array([7.6091, 6.6117, 9.0230]),
        "sym_ops_str": [
            '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2', '-x+1/2,y+1/2,z',
            '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
            'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2', 'x+1/2,y+1/2,z',
            'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
        ],
        "atoms_asymm": [
            {'label': 'Al_A', 'element': 'Al', 'frac_coords': np.array([0.3182, 0.2158, 0.2500])},
            {'label': 'Al_B', 'element': 'Al', 'frac_coords': np.array([0.0000, 0.3662, 0.1030])},
            {'label': 'Al_C', 'element': 'Al', 'frac_coords': np.array([0.1743, 0.0000, 0.0000])},
            {'label': 'Re_A', 'element': 'Re', 'frac_coords': np.array([0.0000, 0.0445, 0.2500])}
        ]
    }

    def parse_sym_op(op_str):
        """Parses a symmetry operation string into a function that can be applied to coordinates."""
        parts = op_str.split(',')
        def apply_op(coords):
            x, y, z = coords[0], coords[1], coords[2]
            local_scope = {'x': x, 'y': y, 'z': z}
            new_x = eval(parts[0], {}, local_scope)
            new_y = eval(parts[1], {}, local_scope)
            new_z = eval(parts[2], {}, local_scope)
            return np.array([new_x, new_y, new_z])
        return apply_op

    cell_lengths = cif_data["cell_lengths"]
    
    def to_cartesian(frac_coords):
        """Converts fractional coordinates to Cartesian coordinates for an orthorhombic cell."""
        return frac_coords * cell_lengths

    # 2. Generate all atoms in the unit cell
    sym_ops = [parse_sym_op(s) for s in cif_data['sym_ops_str']]
    unit_cell_atoms = []
    seen_coords = []
    tolerance = 1e-4

    for atom_info in cif_data['atoms_asymm']:
        for op in sym_ops:
            frac_coords = op(atom_info['frac_coords'])
            frac_coords -= np.floor(frac_coords) # Normalize to [0, 1)
            
            is_duplicate = False
            for seen in seen_coords:
                diff = frac_coords - seen
                diff -= np.round(diff) # Periodic boundary check
                if np.linalg.norm(diff) < tolerance:
                    is_duplicate = True
                    break
            
            if not is_duplicate:
                unit_cell_atoms.append({'element': atom_info['element'], 'frac_coords': frac_coords})
                seen_coords.append(frac_coords)

    # 3. Create a supercell
    supercell_atoms = []
    for atom in unit_cell_atoms:
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    supercell_atoms.append({
                        'element': atom['element'],
                        'frac_coords': atom['frac_coords'] + np.array([i, j, k])
                    })

    # 4. Find coordination for each unique site
    polyhedra_results = []
    
    for central_atom in cif_data['atoms_asymm']:
        central_cart = to_cartesian(central_atom['frac_coords'])
        
        # 4a. Calculate distances
        distances = []
        for neighbor in supercell_atoms:
            neighbor_cart = to_cartesian(neighbor['frac_coords'])
            dist = np.linalg.norm(central_cart - neighbor_cart)
            
            if dist > tolerance: # Exclude the atom itself
                distances.append({'dist': dist, 'element': neighbor['element']})
        
        distances.sort(key=lambda x: x['dist'])
        
        if not distances:
            continue
            
        # 4b. Find coordination shell by largest gap
        first_shell_end_index = len(distances) - 1
        if len(distances) > 1:
            max_search_idx = min(len(distances) - 1, 20)
            gaps = [distances[i+1]['dist'] / distances[i]['dist'] for i in range(max_search_idx)]
            if gaps:
                largest_gap_index = np.argmax(gaps)
                first_shell_end_index = largest_gap_index

        first_shell_neighbors = distances[:first_shell_end_index + 1]
        
        # 4c. Determine formula
        neighbor_elements = [n['element'] for n in first_shell_neighbors]
        counts = Counter(neighbor_elements)
        
        formula_parts = []
        for elem in sorted(counts.keys()):
            count = counts[elem]
            formula_parts.append(f"{elem}{count if count > 1 else ''}")
        
        formula = "".join(formula_parts)
        polyhedra_results.append(f"{central_atom['element']}, {formula}")

    # 5. Format and print the unique results
    unique_results = sorted(list(set(polyhedra_results)))
    print("; ".join(unique_results))

solve()
<<<Al, Al10Re2; Re, Al12>>>