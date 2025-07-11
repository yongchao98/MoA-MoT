import numpy as np
from collections import Counter
import math

def find_coordination_polyhedra():
    """
    This script determines the coordination polyhedra for a crystal structure
    defined by its cell parameters, basis atoms, and symmetry operations.
    """
    # Step 1: Parse CIF Data (hardcoded from the prompt)
    cell_lengths = np.array([7.609100, 6.611700, 9.023000])
    
    # Basis atoms from the asymmetric unit
    basis_atoms = [
        ("Al_A", "Al", np.array([0.3182, 0.2158, 0.2500])),
        ("Al_B", "Al", np.array([0.0000, 0.3662, 0.1030])),
        ("Al_C", "Al", np.array([0.1743, 0.0000, 0.0000])),
        ("Re_A", "Re", np.array([0.0000, 0.0445, 0.2500])),
    ]

    # Symmetry operations for space group Cmcm (No. 63)
    sym_ops_str = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2', '-x+1/2,y+1/2,z',
        '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2', 'x+1/2,y+1/2,z',
        'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]

    # Helper function to parse and apply a symmetry operation string
    def parse_sym_op(op_str):
        ops = op_str.split(',')
        # Using eval is safe here as the input strings are fixed
        fx = lambda x, y, z: eval(ops[0])
        fy = lambda x, y, z: eval(ops[1])
        fz = lambda x, y, z: eval(ops[2])
        return lambda p: np.array([fx(p[0], p[1], p[2]), fy(p[0], p[1], p[2]), fz(p[0], p[1], p[2])])

    sym_ops_func = [parse_sym_op(s) for s in sym_ops_str]

    # Step 2: Generate Full Unit Cell
    all_generated_atoms = []
    for _, symbol, coords in basis_atoms:
        for op_func in sym_ops_func:
            all_generated_atoms.append((symbol, op_func(coords)))

    unique_atoms_in_cell = []
    seen_coords = set()
    for symbol, coords in all_generated_atoms:
        # Normalize coordinates to be within [0, 1) to handle periodic duplicates
        norm_coords = np.mod(coords, 1.0)
        # Round to 4 decimal places to avoid floating point inaccuracies for comparison
        coord_tuple = tuple(np.round(norm_coords, 4))
        if coord_tuple not in seen_coords:
            seen_coords.add(coord_tuple)
            unique_atoms_in_cell.append((symbol, norm_coords))

    # Step 3: Construct a Supercell
    supercell_atoms = []
    for symbol, coords in unique_atoms_in_cell:
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:
                    supercell_atoms.append((symbol, coords + np.array([i, j, k])))

    # Function to calculate Cartesian distance from fractional coordinates
    def get_distance(p1, p2, cell):
        delta_frac = p1 - p2
        delta_cart = delta_frac * cell
        return np.linalg.norm(delta_cart)

    polyhedra = []

    # Steps 4 & 5 & 6: For each basis atom, find neighbors and determine formula
    for _, central_symbol, central_coords in basis_atoms:
        distances = []
        for neighbor_symbol, neighbor_coords in supercell_atoms:
            dist = get_distance(central_coords, neighbor_coords, cell_lengths)
            # Exclude distance to the atom itself (which is 0)
            if dist > 1e-4:
                distances.append((dist, neighbor_symbol))
        
        distances.sort()
        
        if not distances:
            continue
            
        # Identify the first coordination shell using a 15% relative tolerance
        first_dist = distances[0][0]
        cutoff = first_dist * 1.15
        
        neighbors = [symbol for dist, symbol in distances if dist < cutoff]
        
        # Count elements to form the formula
        counts = Counter(neighbors)
        # Sort by element symbol (Al, Re) for a consistent formula string
        formula = "".join(f"{elem}{counts[elem]}" for elem in sorted(counts.keys()))
        
        polyhedra.append(f"{central_symbol}, {formula}")

    # Step 7: Format Output
    # Find the unique polyhedra and sort them for a canonical answer
    unique_polyhedra = sorted(list(set(polyhedra)))
    print("; ".join(unique_polyhedra))

find_coordination_polyhedra()