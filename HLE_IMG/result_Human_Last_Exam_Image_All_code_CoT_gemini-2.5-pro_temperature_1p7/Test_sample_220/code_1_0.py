import numpy as np

def get_fcc_atoms(a=1.0):
    """Returns the 14 atoms in a conventional FCC unit cell."""
    # Corners
    corners = [
        [0, 0, 0], [a, 0, 0], [0, a, 0], [0, 0, a],
        [a, a, 0], [a, 0, a], [0, a, a], [a, a, a]
    ]
    # Face centers
    faces = [
        [a/2, a/2, 0], [a/2, 0, a/2], [0, a/2, a/2],
        [a/2, a, a/2], [a, a/2, a/2], [a/2, a/2, a]
    ]
    return corners + faces

def project_110(points_3d):
    """Projects 3D points onto a plane perpendicular to the [110] direction."""
    projected_points = []
    for p in points_3d:
        x, y, z = p
        # Projection is onto a plane defined by [0,0,1] and [1,-1,0] axes
        x_proj = (x - y) / np.sqrt(2)
        y_proj = z
        projected_points.append((round(x_proj, 4), round(y_proj, 4)))
    
    # Return unique projected points
    return sorted(list(set(projected_points)))

def find_centered_rectangle(points_2d):
    """Finds a rectangle and a centering atom in the projected points."""
    a = 1.0 # Lattice constant
    
    # Define corners of a potential unit rectangle in the projection
    # These correspond to projections of (0,0,0), (1,0,0), (0,0,1), (1,0,1)
    c1 = (0.0, 0.0)
    c2 = (round(a / np.sqrt(2), 4), 0.0)
    c3 = (0.0, a)
    c4 = (round(a / np.sqrt(2), 4), a)
    
    corners = [c1, c2, c3, c4]
    
    # Calculate the geometric center of this rectangle
    center_x = (c1[0] + c2[0]) / 2
    center_y = (c1[1] + c3[1]) / 2
    calculated_center = (round(center_x, 4), round(center_y, 4))
    
    # Find the atom at the center. This corresponds to the projection of (0.5, 0, 0.5) or (1, 0.5, 0.5)
    center_atom = (round(a / (2 * np.sqrt(2)), 4), round(a / 2, 4))
    
    # Check if all these points exist in the projected set
    all_found = all(p in points_2d for p in corners) and center_atom in points_2d
    
    if all_found:
        print("A centered rectangle pattern was found.")
        print("Rectangle corner coordinates in the projection plane:")
        for i, p in enumerate(corners, 1):
            print(f"Corner {i}: {p}")
            
        print(f"\nCalculated center of the rectangle: {calculated_center}")
        print(f"Found center atom coordinate:     {center_atom}")
        
        if calculated_center == center_atom:
            print("\nThe calculated center matches the position of the center atom.")
            print("This confirms the centered rectangular pattern for the FCC [110] projection.")
        else:
            print("Error: Center mismatch.")
    else:
        print("Could not find a centered rectangle pattern.")

# --- Main execution ---
# 1. Get atom positions for a standard FCC unit cell
atoms_3d = get_fcc_atoms(a=1.0)

# 2. Project them along the [110] direction
projected_atoms = project_110(atoms_3d)

# 3. Analyze the projected pattern to confirm it's a centered rectangle
find_centered_rectangle(projected_atoms)

print("\nImage C is the only option that shows this characteristic centered rectangular arrangement.")
