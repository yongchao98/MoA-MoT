import numpy as np

def generate_fcc_projection():
    """
    Generates and analyzes the 2D projection of an FCC lattice
    when viewed along the [110] direction.
    """
    atoms_3d = []
    # Generate points (i, j, k) in a 4x4x4 cube where i+j+k is even.
    # This represents an FCC lattice.
    for i in range(-2, 3):
        for j in range(-2, 3):
            for k in range(-2, 3):
                if (i + j + k) % 2 == 0:
                    atoms_3d.append((i, j, k))

    # Project the 3D atoms onto a 2D plane perpendicular to [110].
    # The new axes are x' = y - x and y' = z.
    projected_atoms = set()
    for x, y, z in atoms_3d:
        xp = y - x
        yp = z
        projected_atoms.add((xp, yp))

    # Find a sample rectangle and its center to demonstrate the pattern
    # The pattern consists of points where (x', y') have the same parity.
    # e.g., (even, even) and (odd, odd)
    
    # Let's find corners of a rectangle. Let a corner be (0,0).
    # Since all coordinates have the same parity, the next horizontal point on the grid is (2,0)
    # The next vertical point is (0,2). So a rectangle has corners (0,0), (2,0), (0,2), (2,2)
    p1 = (0, 0)
    p2 = (2, 0)
    p3 = (0, 2)
    p4 = (2, 2)
    
    # Calculate the center of this rectangle
    center_x = (p1[0] + p4[0]) / 2
    center_y = (p1[1] + p4[1]) / 2
    center_p = (int(center_x), int(center_y))

    print("Analysis of the FCC [110] projection:")
    print("The projection of an FCC lattice along the [110] direction results in a 2D pattern where all coordinates (x', y') have the same parity.")
    print("This creates a 'centered rectangular' lattice pattern.\n")

    print("For example, a unit rectangle in this pattern can be defined by corner points:")
    print(f"Corner 1 = {p1}")
    print(f"Corner 2 = {p2}")
    print(f"Corner 3 = {p3}")
    print(f"Corner 4 = {p4}")
    
    # Check if these points exist in our generated set
    all_found = all(p in projected_atoms for p in [p1,p2,p3,p4,center_p])
    
    if all_found:
        print("\nAll these corner points are present in the projected lattice.")
        print(f"The center of this rectangle is located at the coordinates given by the equation: (({p1[0]} + {p4[0]})/2, ({p1[1]} + {p4[1]})/2)")
        print(f"Which evaluates to: {center_p}")
        print(f"An atom is found at this center point {center_p}, confirming the centered rectangular structure.")
    
    print("\nComparing this theoretical pattern to the given images:")
    print("Image B displays a centered rectangular pattern, consistent with our findings.")
    print("Image A shows a simple rectangular pattern (characteristic of BCC [110]).")
    print("Image C and D show more complex patterns, likely from different structures (like diamond-cubic) or different viewing directions.")
    print("\nTherefore, pattern B indicates a face-centered cubic structure viewed along the [110] direction.")

generate_fcc_projection()
print("\nB.")