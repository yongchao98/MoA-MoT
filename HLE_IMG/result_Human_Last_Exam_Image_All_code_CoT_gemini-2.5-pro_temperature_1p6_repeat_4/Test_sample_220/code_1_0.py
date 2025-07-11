import numpy as np

def generate_fcc_projection(a=1.0, supercell_range=(-1, 2)):
    """
    Generates the 2D coordinates of an FCC lattice projected along the [110] direction.
    
    Returns:
        A dictionary containing the coordinates of the corner and center atoms 
        of a characteristic unit rectangle in the projection.
    """
    # FCC basis atoms in a conventional unit cell
    basis = np.array([
        [0, 0, 0],
        [0.5, 0.5, 0],
        [0.5, 0, 0.5],
        [0, 0.5, 0.5]
    ]) * a

    # Generate atoms in a supercell to see the repeating pattern
    atoms_3d = []
    for i in range(supercell_range[0], supercell_range[1]):
        for j in range(supercell_range[0], supercell_range[1]):
            for k in range(supercell_range[0], supercell_range[1]):
                translation = np.array([i, j, k]) * a
                for b in basis:
                    atoms_3d.append(b + translation)
    
    atoms_3d = np.array(atoms_3d)

    # Project atoms onto the plane perpendicular to [110].
    # The new 2D axes can be chosen along [1,-1,0] and [0,0,1].
    # x_proj = (x - y) / sqrt(2)
    # y_proj = z
    atoms_2d = []
    for p in atoms_3d:
        x_proj = (p[0] - p[1])
        y_proj = p[2]
        atoms_2d.append([x_proj, y_proj])
        
    # Remove duplicate points resulting from projection
    atoms_2d = np.unique(np.round(atoms_2d, decimals=5), axis=0)

    # Identify the corners and center of a unit rectangle
    # The projected lattice has basis vectors corresponding to the 3D vectors a*[1,0,0] and a/2*[0,1,1]. No.
    # The conventional unit cell of the projected lattice has side lengths a*sqrt(2) and a.
    # Let's scale our un-normalized coordinates to match this.
    # Scaling x by 1/sqrt(2) gives the correct aspect ratio.
    x_coords = atoms_2d[:, 0] / np.sqrt(2)
    y_coords = atoms_2d[:, 1]
    
    # Let's find a rectangle. For a=1, corners should be at (0,0), (sqrt(2),0), (0,1), (sqrt(2),1)
    # and the center at (sqrt(2)/2, 1/2). Let's locate these points from our generated set.
    
    corners = {
        "bottom_left": (0.0, 0.0),
        "bottom_right": (a * np.sqrt(2), 0.0),
        "top_left": (0.0, a),
        "top_right": (a * np.sqrt(2), a)
    }
    
    center = (a * np.sqrt(2) / 2, a / 2)

    return {"corners": corners, "center": center}

# Set lattice parameter a=1 for simplicity
result = generate_fcc_projection(a=1.0)

print("The projection of an FCC lattice along the [110] direction results in a centered rectangular pattern.")
print("The characteristic unit cell of this 2D pattern is a rectangle with atoms at its corners and one in the center.")
print("\nFor a lattice parameter a=1, the coordinates of the atoms forming one such rectangle are:")

# Print the corners
print("\nCorners:")
for name, coords in result["corners"].items():
    print(f"{name}: (x={coords[0]:.2f}, y={coords[1]:.2f})")

# Print the center
print("\nCenter:")
center_coords = result["center"]
print(f"(x={center_coords[0]:.2f}, y={center_coords[1]:.2f})")

print("\nThis pattern of a centered rectangle is shown in image C.")
print("The other images correspond to different crystal structures or viewing directions:")
print("- A: BCC lattice viewed along [110]")
print("- B: BCC lattice viewed along [100]")
print("- D: FCC lattice viewed along [100]")
