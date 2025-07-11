import numpy as np
from collections import Counter
import math

def solve():
    # Step 1: Parse Crystal Data from the CIF
    # Cell parameters for ReAl12
    a, b, c = 7.6091, 6.6117, 9.023
    alpha, beta, gamma = 90.0, 90.0, 90.0

    # For an orthorhombic cell, the fractional-to-Cartesian transformation matrix is diagonal
    frac_to_cart_matrix = np.array([[a, 0, 0], [0, b, 0], [0, 0, c]])

    # Asymmetric unit atoms: [label, element_type, [fract_x, fract_y, fract_z]]
    asymmetric_atoms = [
        ["Al_A", "Al", np.array([0.3182, 0.2158, 0.2500])],
        ["Al_B", "Al", np.array([0.0000, 0.3662, 0.1030])],
        ["Al_C", "Al", np.array([0.1743, 0.0000, 0.0000])],
        ["Re_A", "Re", np.array([0.0000, 0.0445, 0.2500])]
    ]

    # Symmetry operations for Cmcm (No. 63)
    sym_ops_str = [
        "-x+1/2,-y+1/2,-z", "-x+1/2,-y+1/2,z+1/2", "-x+1/2,y+1/2,-z+1/2", "-x+1/2,y+1/2,z",
        "-x,-y,-z", "-x,-y,z+1/2", "-x,y,-z+1/2", "-x,y,z",
        "x+1/2,-y+1/2,-z", "x+1/2,-y+1/2,z+1/2", "x+1/2,y+1/2,-z+1/2", "x+1/2,y+1/2,z",
        "x,-y,-z", "x,-y,z+1/2", "x,y,-z+1/2", "x,y,z"
    ]

    # A function to apply a text-based symmetry operation
    def apply_sym_op(coords, op_str):
        x, y, z = coords[0], coords[1], coords[2]
        # Using eval is a concise way to parse the operation strings for this problem
        new_coords = [eval(part) for part in op_str.split(',')]
        return np.array(new_coords)

    # Step 2: Generate the full unit cell
    full_cell_atoms = []
    seen_positions = set()
    for _, element, frac_coords in asymmetric_atoms:
        for op_str in sym_ops_str:
            new_frac = apply_sym_op(frac_coords, op_str)
            # Normalize coordinates to be within the [0, 1) range
            new_frac -= np.floor(new_frac)
            
            # Round to avoid floating point issues and check for uniqueness
            pos_tuple = tuple(np.round(new_frac, 4))
            if pos_tuple not in seen_positions:
                seen_positions.add(pos_tuple)
                full_cell_atoms.append({'element': element, 'frac_coords': new_frac})

    # Step 3, 4, 5: Find coordination polyhedra for each unique atom
    polyhedra_results = set()
    
    for central_label, central_element, central_frac_coords in asymmetric_atoms:
        neighbors = []
        # Search in a 3x3x3 block of cells to find the nearest periodic images
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    for atom in full_cell_atoms:
                        # Calculate the vector from the central atom to the neighbor image
                        vec_frac = atom['frac_coords'] + np.array([dx, dy, dz]) - central_frac_coords
                        
                        # Convert vector to Cartesian coordinates
                        vec_cart = np.dot(frac_to_cart_matrix, vec_frac)
                        
                        # Calculate distance
                        distance = np.linalg.norm(vec_cart)
                        
                        # Store neighbor if it's not the central atom itself (distance > 0.01)
                        if distance > 0.01:
                            neighbors.append({'element': atom['element'], 'distance': distance})
        
        # Sort neighbors by distance
        neighbors.sort(key=lambda x: x['distance'])
        
        if not neighbors:
            continue

        # Identify the first coordination shell by finding a large gap in distances
        first_shell_atoms = []
        if neighbors:
            first_shell_atoms.append(neighbors[0])
            for i in range(1, len(neighbors)):
                # A gap of > 0.5 Angstrom typically indicates the end of the first shell
                if (neighbors[i]['distance'] - neighbors[i-1]['distance']) < 0.5:
                    first_shell_atoms.append(neighbors[i])
                else:
                    break
        
        # Count the elements in the first shell to determine the formula
        shell_elements = [n['element'] for n in first_shell_atoms]
        counts = Counter(shell_elements)
        formula = "".join(f"{el}{counts[el]}" for el in sorted(counts.keys()))
        
        polyhedra_results.add(f"{central_element}, {formula}")

    # Step 6: Format and print the final results
    final_answer = "; ".join(sorted(list(polyhedra_results)))
    print(final_answer)

solve()