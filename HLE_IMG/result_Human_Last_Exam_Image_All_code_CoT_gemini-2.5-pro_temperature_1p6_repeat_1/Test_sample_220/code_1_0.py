import numpy as np

def generate_fcc_projection():
    """
    Generates and projects FCC atomic coordinates for a [110] view.

    The function first defines the basis atoms of an FCC unit cell. It then
    generates atoms in a 2x2x2 supercell. Finally, it projects these atoms
    onto the plane perpendicular to the [110] direction, using the [001]
    and [1,-1,0] directions as the projection axes.

    The projected points are separated into two sets based on their depth to
    illustrate the two interpenetrating rectangular lattices that form the
    final "centered rectangle" pattern.
    """
    # Let the lattice constant a = 1 for simplicity.
    a = 1.0

    # Basis atoms in a conventional FCC unit cell (coordinates as fractions of a)
    basis_atoms = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5]
    ]

    # Generate atoms in a 2x2x2 supercell to see the repeating pattern
    atoms_3d = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for basis_atom in basis_atoms:
                    atoms_3d.append([
                        i + basis_atom[0],
                        j + basis_atom[1],
                        k + basis_atom[2]
                    ])

    # Project atoms onto the plane perpendicular to [110]
    # We use v1 = [1,-1,0] and v2 = [0,0,1] as the axes for the 2D plane.
    # The depth is determined by the [110] direction.
    # projected_x = (x - y) * a
    # projected_y = z * a
    # depth = (x + y) * a
    
    # We use a small tolerance for floating point comparisons
    TOL = 1e-6
    
    layer_A_coords = set()
    layer_B_coords = set()

    for x, y, z in atoms_3d:
        # Check the depth layer. Depth d = x+y.
        # Layers are separated by d=0.5.
        # Layer A: d = integer (0, 1, 2, ...)
        # Layer B: d = half-integer (0.5, 1.5, ...)
        depth = x + y
        
        # We scale the x-projection for a better aspect ratio in the output
        proj_x = round((x - y) * a, 4)
        proj_y = round(z * a, 4)

        # Check if depth is close to an integer
        if abs(depth - round(depth)) < TOL:
            layer_A_coords.add((proj_x, proj_y))
        else: # Depth is half-integer
            layer_B_coords.add((proj_x, proj_y))

    print("The projection of an FCC lattice along the [110] direction forms a centered rectangle.")
    print("This pattern can be seen as two interpenetrating rectangular grids.")
    print("-" * 30)

    print("Projected coordinates for Layer A (forms a rectangular grid):")
    # Sort for consistent output
    sorted_A = sorted(list(layer_A_coords))
    for coord in sorted_A:
        print(f"({coord[0]:>5}, {coord[1]:>5})")
        
    print("\nProjected coordinates for Layer B (forms a centered rectangular grid):")
    sorted_B = sorted(list(layer_B_coords))
    for coord in sorted_B:
        print(f"({coord[0]:>5}, {coord[1]:>5})")
        
    print("-" * 30)
    print("Image B shows this characteristic centered rectangular pattern.")


generate_fcc_projection()