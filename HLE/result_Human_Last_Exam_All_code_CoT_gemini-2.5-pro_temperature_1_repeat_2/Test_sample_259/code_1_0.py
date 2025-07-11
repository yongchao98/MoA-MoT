import numpy as np
from collections import Counter
import math

def solve():
    """
    This function analyzes the provided CIF data for ReAl12 to determine
    the coordination polyhedra for each unique atomic site.
    """
    # Step 1: Define the crystal structure data from the CIF file.
    # This includes cell parameters, symmetry operations, and atomic positions in the asymmetric unit.
    cell_lengths = np.array([7.6091, 6.6117, 9.0230])
    # For orthorhombic cell, angles are 90, so distance calculation is simpler.
    
    # Symmetry operations for space group Cmcm (No. 63)
    sym_ops_str = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2',
        '-x+1/2,y+1/2,z', '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2',
        'x+1/2,y+1/2,z', 'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]

    # Atom sites in the asymmetric unit
    atoms_asu = [
        {'label': 'Al_A', 'symbol': 'Al', 'coords': np.array([0.3182, 0.2158, 0.2500])},
        {'label': 'Al_B', 'symbol': 'Al', 'coords': np.array([0.0000, 0.3662, 0.1030])},
        {'label': 'Al_C', 'symbol': 'Al', 'coords': np.array([0.1743, 0.0000, 0.0000])},
        {'label': 'Re_A', 'symbol': 'Re', 'coords': np.array([0.0000, 0.0445, 0.2500])},
    ]

    # Helper function to apply a symmetry operation string to coordinates
    def apply_sym_op(op_str, coords):
        x, y, z = coords
        # Using a dictionary for eval's local scope is safer than global scope
        local_scope = {'x': x, 'y': y, 'z': z, 'math': math}
        coords_str_list = op_str.split(',')
        new_x = eval(coords_str_list[0], {}, local_scope)
        new_y = eval(coords_str_list[1], {}, local_scope)
        new_z = eval(coords_str_list[2], {}, local_scope)
        return np.array([new_x, new_y, new_z])

    # Step 2: Generate all atom positions in the full unit cell.
    # Apply all symmetry operations to each atom in the asymmetric unit.
    all_atoms_cell = []
    for atom_asu in atoms_asu:
        for op_str in sym_ops_str:
            new_coords = apply_sym_op(op_str, atom_asu['coords'])
            # Normalize coordinates to be within [0, 1) for easier duplicate checking
            new_coords -= np.floor(new_coords)
            all_atoms_cell.append({'symbol': atom_asu['symbol'], 'coords': new_coords})
            
    # Remove duplicate atoms that result from site symmetries
    unique_atoms_cell = []
    for atom in all_atoms_cell:
        is_duplicate = any(
            atom['symbol'] == unique_atom['symbol'] and 
            np.allclose(atom['coords'], unique_atom['coords'], atol=1e-4)
            for unique_atom in unique_atoms_cell
        )
        if not is_duplicate:
            unique_atoms_cell.append(atom)
            
    # Step 3: Create a supercell (3x3x3) to find nearest neighbors correctly,
    # accounting for atoms in adjacent unit cells.
    atoms_supercell = []
    for atom in unique_atoms_cell:
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    atoms_supercell.append({
                        'symbol': atom['symbol'],
                        'coords': atom['coords'] + np.array([i, j, k])
                    })

    # Helper function to calculate distance considering periodic boundary conditions
    def get_distance(p1, p2, cell_lengths):
        # Minimum image convention for periodic boundaries
        delta = p1 - p2
        delta -= np.round(delta)
        # Convert fractional coordinate difference to Cartesian for distance calc
        cart_delta = delta * cell_lengths
        return np.linalg.norm(cart_delta)

    # Step 4: For each unique atom site, find its coordination polyhedron.
    polyhedra_results = []
    for central_atom in atoms_asu:
        # Calculate distances to all other atoms in the supercell
        distances = []
        for neighbor_atom in atoms_supercell:
            # Skip calculating distance to the central atom itself
            if np.allclose(central_atom['coords'], neighbor_atom['coords'], atol=1e-5):
                continue
            
            dist = get_distance(central_atom['coords'], neighbor_atom['coords'], cell_lengths)
            distances.append({'dist': dist, 'symbol': neighbor_atom['symbol']})

        # Sort atoms by distance to find the nearest neighbors
        distances.sort(key=lambda x: x['dist'])
        
        if not distances:
            continue

        # Step 5: Identify the first coordination shell.
        # A common heuristic is to find a large gap between consecutive sorted distances.
        # A gap larger than 0.5 Angstrom is a strong indicator of a new shell.
        cutoff_dist = distances[-1]['dist'] + 1.0 # Default to include all if no gap found
        for i in range(len(distances) - 1):
            gap = distances[i+1]['dist'] - distances[i]['dist']
            if gap > 0.5:
                cutoff_dist = (distances[i]['dist'] + distances[i+1]['dist']) / 2.0
                break
        
        # Collect all neighbors within the first coordination shell
        first_shell_neighbors = [d['symbol'] for d in distances if d['dist'] < cutoff_dist]

        # Step 6: Count the atom types in the shell to create the chemical formula.
        counts = Counter(first_shell_neighbors)
        formula_parts = []
        # Sort by element symbol for a canonical formula (Al, Re)
        for symbol in sorted(counts.keys()):
            count = counts[symbol]
            formula_parts.append(f"{symbol}{count}")
        formula = "".join(formula_parts)
        
        central_symbol = central_atom['symbol']
        polyhedra_results.append(f"{central_symbol}, {formula}")
    
    # Remove duplicate polyhedra definitions and format the final output string.
    unique_results = sorted(list(set(polyhedra_results)))
    print("; ".join(unique_results))

solve()