import numpy as np
from collections import Counter

def get_coordination_polyhedra():
    """
    Analyzes the provided crystal structure data to determine the coordination polyhedra for each unique atomic site.
    """
    # --- Step 1: Parse Crystal Data ---
    # Lattice parameters for ReAl12
    a, b, c = 7.6091, 6.6117, 9.0230
    
    # Transformation matrix for an orthorhombic cell (converts fractional to Cartesian coordinates)
    frac_to_cart_matrix = np.array([
        [a, 0, 0],
        [0, b, 0],
        [0, 0, c]
    ])

    # Symmetry operations for space group Cmcm (No. 63)
    symmetry_operations = [
        '-x+1/2,-y+1/2,-z', '-x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2', '-x+1/2,y+1/2,z',
        '-x,-y,-z', '-x,-y,z+1/2', '-x,y,-z+1/2', '-x,y,z',
        'x+1/2,-y+1/2,-z', 'x+1/2,-y+1/2,z+1/2', 'x+1/2,y+1/2,-z+1/2', 'x+1/2,y+1/2,z',
        'x,-y,-z', 'x,-y,z+1/2', 'x,y,-z+1/2', 'x,y,z'
    ]

    # Unique atom sites in the asymmetric unit (Type, Fractional Coords)
    asymmetric_unit_atoms = [
        ("Re", np.array([0.0000, 0.0445, 0.2500])),
        ("Al", np.array([0.3182, 0.2158, 0.2500])),
        ("Al", np.array([0.0000, 0.3662, 0.1030])),
        ("Al", np.array([0.1743, 0.0000, 0.0000]))
    ]

    # Helper function to apply a symmetry operation to a coordinate
    def apply_symop(coords, op_str):
        x, y, z = coords
        # Using eval is safe here as the operation strings are fixed and well-defined
        op_x, op_y, op_z = op_str.split(',')
        new_x = eval(op_x, {'x': x, 'y': y, 'z': z})
        new_y = eval(op_y, {'x': x, 'y': y, 'z': z})
        new_z = eval(op_z, {'x': x, 'y': y, 'z': z})
        return np.array([new_x, new_y, new_z])

    # --- Step 2: Generate Full Unit Cell ---
    unit_cell_atoms = []
    for atom_type, frac_coords in asymmetric_unit_atoms:
        for op in symmetry_operations:
            new_coords = apply_symop(frac_coords, op)
            # Normalize coordinates into the unit cell [0, 1) range
            new_coords = np.mod(new_coords, 1.0)
            
            # Avoid adding duplicate atoms that result from site symmetry
            is_duplicate = any(np.allclose(new_coords, existing_coords, atol=1e-4) 
                               for _, existing_coords in unit_cell_atoms)
            if not is_duplicate:
                unit_cell_atoms.append((atom_type, new_coords))

    # --- Step 3: Construct a Supercell ---
    # A 3x3x3 block of unit cells ensures all nearest neighbors are found
    supercell_atoms = []
    for atom_type, frac_coords in unit_cell_atoms:
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    supercell_atoms.append((atom_type, frac_coords + np.array([i, j, k])))

    # --- Step 4, 5, 6: Find Distances, Shells, and Formulas ---
    final_polyhedra = []
    for central_atom_type, central_frac_coords in asymmetric_unit_atoms:
        
        distances = []
        for neighbor_type, neighbor_frac_coords in supercell_atoms:
            # Get vector from central atom to neighbor in fractional space
            frac_diff = neighbor_frac_coords - central_frac_coords
            
            # Apply minimum image convention to find the closest periodic image
            frac_diff -= np.round(frac_diff)
            
            # Convert vector to Cartesian coordinates and find its length (the distance)
            cart_diff = frac_to_cart_matrix @ frac_diff
            distance = np.linalg.norm(cart_diff)

            # Exclude the central atom itself (distance is zero)
            if distance > 1e-4:
                distances.append((distance, neighbor_type))
        
        distances.sort()

        # Identify the first coordination shell by finding the first major gap in distances
        if not distances:
            continue
        
        first_shell = []
        last_distance_in_shell = distances[0][0]
        first_shell.append(distances[0])

        # A gap is defined as an increase in distance of >25% from the previous neighbor
        for dist, ntype in distances[1:]:
            if (dist - last_distance_in_shell) > (last_distance_in_shell * 0.25):
                break
            first_shell.append((dist, ntype))
            last_distance_in_shell = dist
        
        # Count atom types in the shell to create the formula
        neighbor_counts = Counter(neighbor[1] for neighbor in first_shell)
        
        # Sort elements alphabetically (Al, Re) for a consistent formula string
        formula = "".join(f"{elem}{count}" for elem, count in sorted(neighbor_counts.items()))
        
        final_polyhedra.append(f"{central_atom_type}, {formula}")

    # --- Step 7: Format and Print the Final Output ---
    # Print each unique polyhedron, sorted for a consistent and predictable result
    print("; ".join(sorted(final_polyhedra)))

if __name__ == '__main__':
    get_coordination_polyhedra()