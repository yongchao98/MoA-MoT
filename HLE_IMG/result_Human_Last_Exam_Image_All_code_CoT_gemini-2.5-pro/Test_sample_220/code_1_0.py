import numpy as np

def project_fcc_110():
    """
    Calculates the projection of FCC atoms along the [110] direction
    and identifies the resulting pattern.
    """
    print("Step 1: Define key atom positions in an FCC lattice (unit cell side a=1).")
    # Corner atoms and face-center atoms that define the projected pattern
    # A corner atom at the origin
    corner_atom_1 = np.array([0, 0, 0])
    # A corner atom that forms one side of the projected rectangle
    corner_atom_2 = np.array([1, 0, 0])
    # A corner atom that forms the other side of the projected rectangle
    corner_atom_3 = np.array([0, 0, 1])
    # The fourth corner of the projected rectangle
    corner_atom_4 = np.array([1, 0, 1])
    # A face-center atom that projects to the center of the rectangle
    center_atom = np.array([0.5, 0, 0.5])
    # Note: We use un-normalized projection vectors for simplicity. The resulting
    #       pattern shape is what matters. The vector [1,-1,0] has length sqrt(2).
    #       We will scale the x-coordinates by 1/sqrt(2) for a proper aspect ratio.
    scaling_factor = 1 / np.sqrt(2)


    print("\nStep 2: Define the projection axes for the [110] view.")
    # The new 2D axes are perpendicular to the view direction [1, 1, 0]
    proj_axis_x = np.array([1, -1, 0])
    proj_axis_y = np.array([0, 0, 1])
    print(f"x' axis is along {proj_axis_x}, y' axis is along {proj_axis_y}.")

    print("\nStep 3: Calculate the projected 2D coordinates of these key atoms.")
    atoms = {
        "Rectangle Corner 1": corner_atom_1,
        "Rectangle Corner 2": corner_atom_2,
        "Rectangle Corner 3": corner_atom_3,
        "Rectangle Corner 4": corner_atom_4,
        "Center Atom": center_atom
    }
    
    print("The projected coordinates (x', y') form a centered rectangle:")
    for name, pos_3d in atoms.items():
        # Project the 3D position onto the 2D plane
        pos_2d_x = np.dot(pos_3d, proj_axis_x) * scaling_factor
        pos_2d_y = np.dot(pos_3d, proj_axis_y)
        print(f"- {name}: Original 3D position = {pos_3d}")
        print(f"  Projected 2D position (x', y') = ({pos_2d_x:.2f}, {pos_2d_y:.2f})")

    print("\nStep 4: Analyze the pattern and compare with the given images.")
    print("The calculated pattern is a 'centered rectangle'. The atoms form a rectangular grid,")
    print("with another set of atoms located in the center of each rectangle.")
    print("\nVisual comparison:")
    print(" - Image A shows a simple rectangular lattice.")
    print(" - Image B shows a staggered row pattern.")
    print(" - Image C shows a distorted hexagonal-like pattern.")
    print(" - Image D shows a clear centered rectangular lattice, matching our calculated pattern.")
    print("\nConclusion: Image D correctly represents the FCC structure viewed along the [110] direction.")

project_fcc_110()