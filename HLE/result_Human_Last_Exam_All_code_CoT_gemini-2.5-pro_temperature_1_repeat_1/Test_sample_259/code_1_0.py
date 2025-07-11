import numpy as np
from collections import Counter
import re

def solve_coordination():
    """
    Analyzes a crystal structure from a CIF string to find the coordination polyhedra for each unique atom.
    """
    cif_string = """
_cell_length_a    7.609100
_cell_length_b    6.611700
_cell_length_c    9.023000
_cell_angle_alpha 90.000000
_cell_angle_beta  90.000000
_cell_angle_gamma 90.000000
_symmetry_Int_Tables_number    63
_symmetry_space_group_name_H-M Cmcm

loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
 1 -x+1/2,-y+1/2,-z
 2 -x+1/2,-y+1/2,z+1/2
 3 -x+1/2,y+1/2,-z+1/2
 4 -x+1/2,y+1/2,z
 5 -x,-y,-z
 6 -x,-y,z+1/2
 7 -x,y,-z+1/2
 8 -x,y,z
 9 x+1/2,-y+1/2,-z
 10 x+1/2,-y+1/2,z+1/2
 11 x+1/2,y+1/2,-z+1/2
 12 x+1/2,y+1/2,z
 13 x,-y,-z
 14 x,-y,z+1/2
 15 x,y,-z+1/2
 16 x,y,z

loop_
 _atom_site_label
 _atom_site_type_symbol
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 Al_A      Al     0.3182  0.2158  0.2500
 Al_B      Al     0.0000  0.3662  0.1030
 Al_C      Al     0.1743  0.0000  0.0000
 Re_A      Re     0.0000  0.0445  0.2500
"""

    def parse_cif(cif_text):
        data = {
            'cell_lengths': [],
            'sym_ops': [],
            'atoms_au': []
        }
        lines = cif_text.strip().split('\n')
        in_sym_loop = False
        in_atom_loop = False

        for line in lines:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
            
            parts = re.split(r'\s+', line)
            
            if parts[0] == '_cell_length_a': data['cell_lengths'].append(float(parts[1]))
            elif parts[0] == '_cell_length_b': data['cell_lengths'].append(float(parts[1]))
            elif parts[0] == '_cell_length_c': data['cell_lengths'].append(float(parts[1]))
            elif line.startswith('loop_'):
                in_sym_loop, in_atom_loop = False, False
            elif line == '_symmetry_equiv_pos_as_xyz':
                in_sym_loop, in_atom_loop = True, False
            elif line == '_atom_site_label':
                in_atom_loop, in_sym_loop = True, False
            elif in_sym_loop and not parts[0].startswith('_'):
                op_string = line.split(" ", 1)[-1].strip().replace("'", "")
                data['sym_ops'].append(op_string)
            elif in_atom_loop and not parts[0].startswith('_'):
                data['atoms_au'].append({
                    'label': parts[0],
                    'symbol': parts[1],
                    'coords': np.array([float(parts[2]), float(parts[3]), float(parts[4])])
                })
        return data

    def apply_sym_op(coords, op_string):
        x, y, z = coords
        op_parts = op_string.lower().split(',')
        new_x = eval(op_parts[0])
        new_y = eval(op_parts[1])
        new_z = eval(op_parts[2])
        return np.array([new_x, new_y, new_z])

    cif_data = parse_cif(cif_string)
    cell_lengths = np.array(cif_data['cell_lengths'])

    # Generate all atoms in the primary unit cell
    full_cell_atoms = []
    unique_coords_set = set()
    for atom_au in cif_data['atoms_au']:
        for op in cif_data['sym_ops']:
            new_coords_raw = apply_sym_op(atom_au['coords'], op)
            new_coords_mod = np.mod(new_coords_raw, 1.0)
            coord_tuple = tuple(np.round(c, 5) for c in new_coords_mod)
            if coord_tuple not in unique_coords_set:
                unique_coords_set.add(coord_tuple)
                full_cell_atoms.append({'symbol': atom_au['symbol'], 'coords': new_coords_mod})

    # Generate a 3x3x3 supercell for neighbor finding
    supercell_atoms = []
    for atom in full_cell_atoms:
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    translation = np.array([i, j, k])
                    supercell_atoms.append({'symbol': atom['symbol'], 'coords': atom['coords'] + translation})

    results = []
    for central_atom in cif_data['atoms_au']:
        central_coords = central_atom['coords']
        
        distances = []
        for neighbor in supercell_atoms:
            delta_frac = neighbor['coords'] - central_coords
            dist = np.linalg.norm(delta_frac * cell_lengths)

            if dist > 0.01: # Avoid the central atom itself
                distances.append({'dist': dist, 'symbol': neighbor['symbol']})

        distances.sort(key=lambda x: x['dist'])

        if not distances:
            continue

        # Find the largest gap in distances to identify the first coordination shell
        max_gap = 0
        cutoff_index = len(distances) - 1 
        if len(distances) > 1:
            for i in range(len(distances) - 1):
                gap = distances[i+1]['dist'] - distances[i]['dist']
                if gap > 0.1 and gap > max_gap:
                    max_gap = gap
                    cutoff_index = i

        first_shell = distances[:cutoff_index + 1]
        
        # Count elements in the first shell and format the formula
        neighbor_symbols = [d['symbol'] for d in first_shell]
        counts = Counter(neighbor_symbols)
        formula = "".join(f"{elem}{count}" for elem, count in sorted(counts.items()))
        
        results.append(f"{central_atom['symbol']}, {formula}")

    print("; ".join(results))

solve_coordination()