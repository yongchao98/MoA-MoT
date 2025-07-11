import numpy as np

def generate_fcc_projection():
    """
    Generates and prints the 2D coordinates of FCC atoms
    projected along the [110] direction.
    """
    # Let the lattice constant be 'a'
    a = 1.0

    # Define the 14 atoms in a conventional FCC unit cell
    atoms_3d = [
        # 8 corners
        (0, 0, 0), (a, 0, 0), (0, a, 0), (0, 0, a),
        (a, a, 0), (a, 0, a), (0, a, a), (a, a, a),
        # 6 face centers
        (a/2, a/2, 0), (a/2, 0, a/2), (0, a/2, a/2),
        (a, a/2, a/2), (a/2, a, a/2), (a/2, a/2, a)
    ]

    # The [110] viewing direction means we project onto a plane
    # perpendicular to the vector v = [1, 1, 0].
    # We can use the orthogonal vectors u1 = [1, -1, 0] and u2 = [0, 0, 1]
    # as the axes for our 2D projection.
    # A 3D point P=(x,y,z) is projected to (P . u1, P . u2) after normalization.
    # For simplicity, we use (x-y, z) as our projected coordinates.

    projected_atoms = set()
    for x, y, z in atoms_3d:
        # Project the 3D point (x,y,z) to a 2D point (x_proj, z_proj)
        x_proj = round(x - y, 4)
        z_proj = round(z, 4)
        projected_atoms.add((x_proj, z_proj))

    print("Theoretical 2D coordinates for FCC viewed along [110] (lattice constant a=1):")
    # Sort the coordinates for a clear, predictable output
    sorted_coords = sorted(list(projected_atoms))
    for coord in sorted_coords:
        print(f"({coord[0]:.2f}, {coord[1]:.2f})")

    print("\nAnalysis of the pattern:")
    print("The projected points form a centered rectangular lattice.")
    print("For example, consider the rectangle with corners at (-1.00, 0.00), (1.00, 0.00), (-1.00, 1.00), and (1.00, 1.00).")
    print("The pattern has atoms at the corners of this rectangle (e.g. at (1.00, 0.00) and (1.00, 1.00))")
    print("and also at the centers of the short edges (e.g. at (0.00, 0.00) and (0.00, 1.00)).")
    print("Crucially, it also features atoms at the center of the rectangle faces, like at (-0.50, 0.50) and (0.50, 0.50).")
    print("This overall pattern is a centered rectangle.")
    print("\nComparing this theoretical pattern to the images, Image B is the schematic that represents this centered rectangular arrangement for an FCC structure.")


generate_fcc_projection()