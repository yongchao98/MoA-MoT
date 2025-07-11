import math

def solve_candy_packing():
    """
    Solves the candy packing problem by constructing and validating a
    dense, staggered-layer packing arrangement.
    """
    # Box and candy parameters from the problem description
    x_range = (4, 20)
    y_range = (4, 20)
    z_range = (4, 18)
    # The non-overlapping constraint is (xi-xj)^2 + ... >= (4+4)^2 = 64
    min_dist_sq = (4 + 4)**2

    centers = []

    # Constructing the proposed packing with 22 candies in 3 layers
    
    # Layer 1 (z=4): 3x3 grid of 9 candies
    # x,y coordinates: {4, 12, 20}
    for i in range(3):
        for j in range(3):
            centers.append((4 + i * 8, 4 + j * 8, 4))
            
    # Layer 2 (z=11): 2x2 grid of 4 candies, shifted to fit in the hollows
    # x,y coordinates: {8, 16}
    for i in range(2):
        for j in range(2):
            centers.append((8 + i * 8, 8 + j * 8, 11))
            
    # Layer 3 (z=18): 3x3 grid of 9 candies, aligned with Layer 1
    # x,y coordinates: {4, 12, 20}
    for i in range(3):
        for j in range(3):
            centers.append((4 + i * 8, 4 + j * 8, 18))
            
    # Verification Step 1: Check if all centers are within the defined boundaries
    for x, y, z in centers:
        if not (x_range[0] <= x <= x_range[1] and
                y_range[0] <= y <= y_range[1] and
                z_range[0] <= z <= z_range[1]):
            print(f"Error: Center ({x},{y},{z}) is out of bounds.")
            return

    # Verification Step 2: Check non-overlapping constraint for all pairs
    is_valid = True
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            p1 = centers[i]
            p2 = centers[j]
            dist_sq = (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2
            if dist_sq < min_dist_sq:
                print(f"Error: Overlap between {p1} and {p2}. Dist^2 = {dist_sq}")
                is_valid = False
                break
        if not is_valid:
            break

    if is_valid:
        print("The constructed packing is valid.")
        
        # Displaying an example calculation as requested
        p_first = centers[0]         # (4, 4, 4) from layer 1
        p_middle = centers[9]        # (8, 8, 11) from layer 2
        
        x1, y1, z1 = p_first
        x2, y2, z2 = p_middle
        
        dx2 = (x1 - x2)**2
        dy2 = (y1 - y2)**2
        dz2 = (z1 - z2)**2
        total_dist_sq = dx2 + dy2 + dz2
        
        print("\n--- Example Non-Overlapping Check ---")
        print(f"Checking distance between candy at ({x1}, {y1}, {z1}) and ({x2}, {y2}, {z2}):")
        print(f"Equation: ({x1}-{x2})^2 + ({y1}-{y2})^2 + ({z1}-{z2})^2 >= (4+4)^2")
        print(f"Calculation: {dx2} + {dy2} + {dz2} = {total_dist_sq}")
        print(f"Result: {total_dist_sq} is greater than or equal to {min_dist_sq}, so they do not overlap.")
        
        print("\n--- Final Answer ---")
        print(f"The maximized value for n is: {len(centers)}")

solve_candy_packing()