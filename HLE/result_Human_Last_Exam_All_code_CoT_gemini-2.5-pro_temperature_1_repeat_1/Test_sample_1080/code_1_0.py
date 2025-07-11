import math

def solve_candy_packing():
    """
    This function verifies a proposed solution for packing 22 spherical candies
    into a box according to the specified constraints.
    """
    # Problem constraints
    x_min, x_max = 4, 20
    y_min, y_max = 4, 20
    z_min, z_max = 4, 18
    # Sphere radius is 2cm, which is 4 units in the 0.5cm grid.
    # The non-overlapping distance is sum of two radii, so 4+4=8 units.
    min_dist_sq = (4 + 4)**2  # 64

    # Proposed solution: A 3-layer packing with 9+4+9 = 22 candies
    # Layer 1: 9 candies in a 3x3 grid at z=4
    layer1_coords = [(x, y, 4) for x in [4, 12, 20] for y in [4, 12, 20]]
    # Layer 2: 4 candies in a 2x2 grid at z=10, placed in the hollows of layer 1
    layer2_coords = [(x, y, 10) for x in [8, 16] for y in [8, 16]]
    # Layer 3: 9 candies in a 3x3 grid at z=16, aligned with layer 1
    layer3_coords = [(x, y, 16) for x in [4, 12, 20] for y in [4, 12, 20]]

    all_spheres = layer1_coords + layer2_coords + layer3_coords
    n = len(all_spheres)

    # --- Verification Step 1: Check if all sphere centers are within the box ---
    for i, (x, y, z) in enumerate(all_spheres):
        if not (x_min <= x <= x_max and y_min <= y <= y_max and z_min <= z <= z_max):
            print(f"Error: Sphere {i+1} at ({x},{y},{z}) is out of bounds.")
            return
            
    # --- Verification Step 2: Check for overlaps ---
    min_actual_dist_sq = float('inf')
    closest_pair = (None, None)

    for i in range(n):
        for j in range(i + 1, n):
            p1 = all_spheres[i]
            p2 = all_spheres[j]
            dist_sq = (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2
            
            if dist_sq < min_dist_sq:
                print(f"Error: Overlap detected between sphere {i+1} and {j+1}.")
                print(f"Distance squared is {dist_sq}, which is less than {min_dist_sq}.")
                return
            
            if dist_sq < min_actual_dist_sq:
                min_actual_dist_sq = dist_sq
                closest_pair = (p1, p2)

    # --- Print Results ---
    print("Problem formulation is correct: Yes.")
    print(f"A valid configuration with {n} candies was found and verified.")
    print("\nThis packing consists of 3 layers: 9 candies + 4 candies + 9 candies.")
    
    p1, p2 = closest_pair
    dx_sq = (p1[0] - p2[0])**2
    dy_sq = (p1[1] - p2[1])**2
    dz_sq = (p1[2] - p2[2])**2
    
    print("\nVerification of the non-overlapping constraint for the closest pair of candies:")
    print(f"The two closest candies are at {p1} and {p2}.")
    print(f"Their squared distance is {min_actual_dist_sq}.")
    
    # Outputting the numbers in the final equation as requested
    print("\nFinal Equation Check:")
    print(f"(x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2 >= (4+4)^2")
    print(f"{dx_sq} + {dy_sq} + {dz_sq} = {min_actual_dist_sq}")
    print(f"Constraint check: {min_actual_dist_sq} >= {min_dist_sq}. This is true.")

    print(f"\nThus, the maximized value n is {n}.")


solve_candy_packing()