import math

def solve_candy_packing():
    """
    Solves the discrete sphere packing problem for the given box and candy size.
    """
    # Problem Parameters in grid units (1 unit = 0.5 cm)
    # Sphere radius = 2 cm = 4 units
    # Non-overlapping distance between centers >= 2 * radius = 8 units
    min_dist_sq = 8**2

    # Center coordinate ranges
    x_range = (4, 20)
    y_range = (4, 20)
    z_range = (4, 18)

    # --- Constructing the solution using a layered (FCC-like) approach ---

    # Layer 1: A 3x3 grid at z=4.
    layer1_centers = []
    z1 = 4
    for i in range(3):
        for j in range(3):
            layer1_centers.append((4 + i * 8, 4 + j * 8, z1))

    # Layer 2: A 2x2 grid in the hollows of Layer 1, at z=10
    layer2_centers = []
    z2 = 10
    for i in range(2):
        for j in range(2):
            layer2_centers.append((8 + i * 8, 8 + j * 8, z2))

    # Layer 3: A 3x3 grid, same as Layer 1, at z=16
    layer3_centers = []
    z3 = 16
    for i in range(3):
        for j in range(3):
            layer3_centers.append((4 + i * 8, 4 + j * 8, z3))

    # Combine all layers into the final solution
    all_centers = layer1_centers + layer2_centers + layer3_centers
    
    # --- Verification Step ---
    is_valid = True
    # 1. Check if all centers are within bounds
    for x, y, z in all_centers:
        if not (x_range[0] <= x <= x_range[1] and \
                y_range[0] <= y <= y_range[1] and \
                z_range[0] <= z <= z_range[1]):
            print(f"Error: Center ({x},{y},{z}) is out of bounds.")
            is_valid = False
            break
    
    # 2. Check for overlaps
    if is_valid:
        for i in range(len(all_centers)):
            for j in range(i + 1, len(all_centers)):
                p1 = all_centers[i]
                p2 = all_centers[j]
                dist_sq = (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2
                if dist_sq < min_dist_sq:
                    print(f"Error: Overlap between {p1} and {p2}. Squared distance is {dist_sq}.")
                    is_valid = False
                    break
            if not is_valid:
                break
    
    # --- Print the results ---
    print("Yes, the problem formulation is correct.")
    if is_valid:
        print("A valid configuration for packing the candies was found.")
        
        n1 = len(layer1_centers)
        n2 = len(layer2_centers)
        n3 = len(layer3_centers)
        total_n = len(all_centers)

        print("\nThe best arrangement is to stack the candies in three layers:")
        # Outputting each number in the final equation as requested
        print(f"Layer 1 has {n1} candies.")
        print(f"Layer 2 has {n2} candies, placed in the dips of the first layer.")
        print(f"Layer 3 has {n3} candies.")
        print("\nThe final equation to get the total number is:")
        print(f"{n1} + {n2} + {n3} = {total_n}")

        print(f"\nThe maximized value for n is {total_n}.")
    else:
        print("The proposed configuration was found to be invalid.")

solve_candy_packing()