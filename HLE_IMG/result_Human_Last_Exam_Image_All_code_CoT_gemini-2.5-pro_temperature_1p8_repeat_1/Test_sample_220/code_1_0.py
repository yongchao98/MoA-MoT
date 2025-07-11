import numpy as np

def get_fcc_projection_pattern():
    """
    Calculates the projected coordinates of FCC atoms viewed along the [110] direction
    to demonstrate the resulting centered rectangular pattern.
    """
    # Define atoms in an FCC conventional unit cell (lattice constant a=1)
    # These are 8 corners and 6 face centers.
    atoms_3d = [
        # Corners
        (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
        (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1),
        # Face Centers
        (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5),
        (1, 0.5, 0.5), (0.5, 1, 0.5), (0.5, 0.5, 1)
    ]

    # For a [110] view, the 2D projected coordinates are based on (x-y) and z.
    # We select a set of atoms that form one rectangular unit in the projection.
    
    # These four atoms form the corners of a rectangle in the projection.
    rect_corners_3d = [(0, 0, 0), (1, 0, 0), (0, 0, 1), (1, 0, 1)]
    
    # This atom projects to the center of the above rectangle.
    center_atom_3d = (0.5, 0, 0.5)

    print("The projection of an FCC lattice along [110] results in a centered rectangular pattern.")
    print("Here are the 3D coordinates of atoms forming one such unit, and their 2D projection.")
    print("The 2D projection is calculated as (x-y, z).\n")

    print("Rectangle Corner Atoms:")
    for atom in rect_corners_3d:
        x, y, z = atom
        proj_x = x - y
        proj_y = z
        print(f"  Atom at ({x}, {y}, {z}) projects to ({proj_x:.1f}, {proj_y:.1f})")

    print("\nCentering Atom:")
    x, y, z = center_atom_3d
    proj_x = x - y
    proj_y = z
    print(f"  Atom at ({x:.1f}, {y:.1f}, {z:.1f}) projects to ({proj_x:.1f}, {proj_y:.1f})")

    print("\nThis set of projected points forms a rectangle with corners at (0,0), (1,0), (0,1), (1,1) and a centering point at (0.5, 0.5).")
    print("This confirms the 'centered rectangular' nature of the pattern.")
    print("\nVisually, image B is the only one showing a centered pattern.")
    print("Therefore, B is the correct answer.")

get_fcc_projection_pattern()