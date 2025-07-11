import numpy as np
from collections import Counter
import math

def get_coordination_polyhedra():
    """
    This script determines the coordination polyhedra for the ReAl12 crystal structure.
    It parses embedded CIF data, generates all atoms in the unit cell, calculates
    interatomic distances, identifies the first coordination shell using a distance
    gap algorithm, and reports the resulting polyhedra.
    """
    # 1. Parse CIF Data
    # Lattice parameters for orthorhombic cell
    a, b, c = 7.6091, 6.6117, 9.0230
    lattice_matrix = np.array([[a, 0, 0], [0, b, 0], [0, 0, c]])

    # Basis atoms from CIF file
    basis_atoms = {
        'Al_A': ('Al', np.array([0.3182, 0.2158, 0.2500])),
        'Al_B': ('Al', np.array([0.0000, 0.3662, 0.1030])),
        'Al_C': ('Al', np.array([0.1743, 0.0000, 0.0000])),
        'Re_A': ('Re', np.array([0.0000, 0.0445, 0.2500]))
    }

    # Symmetry operations for space group Cmcm (No. 63)
    symops_xyz = [
        '-x+1/2,-y+1/2,-z',  '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2', '-x+1/2,y+1/2,z',
        '-x,-y,-z',         '-x,-y,z+1/2',         '-x,y,-z+1/2',         '-x,y,z',
        'x+1/2,-y+1/2,-z',  'x+1/2,-y+1/2,z+1/2',  'x+1/2,y+1/2,-z+1/2',  'x+1/2,y+1/2,z',
        'x,-y,-z',          'x,-y,z+1/2',          'x,y,-z+1/2',          'x,y,z'
    ]
    
    op_funcs = [eval(f"lambda x, y, z: ({op_str})") for op_str in symops_xyz]

    # 2. Generate the Full Unit Cell
    all_atoms = []
    for site_label, (atom_type, pos) in basis_atoms.items():
        for op in op_funcs:
            new_pos = np.array(op(pos[0], pos[1], pos[2]))
            # Keep coordinates in [0, 1) range
            all_atoms.append({'type': atom_type, 'pos': new_pos - np.floor(new_pos)})

    # Remove duplicate atoms that result from applying symmetry to special positions
    unique_atoms_in_cell = []
    for atom in all_atoms:
        is_duplicate = any(
            atom['type'] == ua['type'] and np.allclose(atom['pos'], ua['pos'], atol=1e-4)
            for ua in unique_atoms_in_cell
        )
        if not is_duplicate:
            unique_atoms_in_cell.append(atom)
    
    results = []
    
    # 3. Find Nearest Neighbors for each unique basis atom
    for center_label, (center_type, center_pos) in basis_atoms.items():
        distances = []
        # Consider a 3x3x3 supercell for periodic boundary conditions
        for neighbor_atom in unique_atoms_in_cell:
            for i in range(-1, 2):
                for j in range(-1, 2):
                    for k in range(-1, 2):
                        # Calculate vector from center to neighbor image
                        delta_frac = neighbor_atom['pos'] + np.array([i, j, k]) - center_pos
                        
                        # Convert fractional vector to Cartesian and get distance
                        dist = np.linalg.norm(np.dot(delta_frac, lattice_matrix))

                        # Exclude distance to the central atom itself (dist ~ 0)
                        if dist > 1e-3:
                            distances.append({'dist': dist, 'type': neighbor_atom['type']})
        
        # Sort neighbors by distance
        distances.sort(key=lambda x: x['dist'])

        # 4. Identify the First Coordination Shell using a distance gap
        if not distances:
            continue
            
        first_shell = []
        last_dist = distances[0]['dist']
        first_shell.append(distances[0])

        for i in range(1, len(distances)):
            dist = distances[i]['dist']
            # If the gap to the previous atom is large, the shell ends
            if dist - last_dist > 0.4:
                break
            first_shell.append(distances[i])
            last_dist = dist
            
        # 5. Formulate and Report the Polyhedra
        neighbor_counts = Counter(atom['type'] for atom in first_shell)
        
        formula = ""
        for elem in sorted(neighbor_counts.keys()):
            count = neighbor_counts[elem]
            formula += elem
            if count > 1:
                formula += str(count)
        
        results.append({
            'label': center_label,
            'type': center_type,
            'formula': f"{center_type}{formula}"
        })

    # Group results by atom type to decide on using labels
    grouped_by_type = {}
    for res in results:
        grouped_by_type.setdefault(res['type'], []).append(res)
    
    final_output = []
    sorted_types = sorted(grouped_by_type.keys())
    
    for atom_type in sorted_types:
        records = grouped_by_type[atom_type]
        if len(records) == 1:
            # If only one unique site for this element, don't use site label
            final_output.append(f"{records[0]['type']}, {records[0]['formula']}")
        else:
            # If multiple sites, distinguish them with labels
            for record in sorted(records, key=lambda x: x['label']):
                 final_output.append(f"{record['label']}, {record['formula']}")

    print("; ".join(final_output))

get_coordination_polyhedra()