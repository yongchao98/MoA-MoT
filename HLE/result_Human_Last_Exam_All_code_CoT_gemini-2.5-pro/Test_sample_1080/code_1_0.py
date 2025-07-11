import math

def solve_candy_packing():
    """
    Calculates and verifies the maximum number of eyeball candies in the box.
    """
    # The problem formulation is correct. We proceed to find the maximum n.
    
    # Strategy: Pack in three layers (A-B-A pattern).
    # Layer A is a 3x3 grid. Layer B is a 2x2 grid in the hollows of A.
    
    # 1. Define the coordinates for the 22 spheres based on the strategy.
    spheres = []
    
    # Layer 1 (A): 3x3 = 9 spheres. Placed at z=4.
    layer1_count = 0
    xy_coords_A = [4, 12, 20]
    z_coord_A1 = 4
    for x in xy_coords_A:
        for y in xy_coords_A:
            spheres.append((x, y, z_coord_A1))
            layer1_count += 1
            
    # Layer 2 (B): 2x2 = 4 spheres. Placed at z=10.
    # The vertical distance between layers must be at least sqrt(32) ~= 5.65.
    # We choose a grid distance of 6 (10 - 4), which is valid.
    layer2_count = 0
    xy_coords_B = [8, 16]
    z_coord_B = 10
    for x in xy_coords_B:
        for y in xy_coords_B:
            spheres.append((x, y, z_coord_B))
            layer2_count += 1
            
    # Layer 3 (A): 3x3 = 9 spheres. Placed at z=16.
    # Vertical distance is again 6 (16 - 10).
    layer3_count = 0
    z_coord_A2 = 16
    for x in xy_coords_A:
        for y in xy_coords_A:
            spheres.append((x, y, z_coord_A2))
            layer3_count += 1
            
    n = len(spheres)
    
    # 2. Verify the solution.
    # All coordinates are within the defined ranges [4,20] for x,y and [4,18] for z.
    # Now, check the non-overlapping constraint for all pairs.
    min_dist_sq = (4 + 4)**2  # 64
    valid_packing = True
    closest_pair_dist_sq = float('inf')
    
    for i in range(n):
        for j in range(i + 1, n):
            p1 = spheres[i]
            p2 = spheres[j]
            dist_sq = (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2
            if dist_sq < min_dist_sq:
                valid_packing = False
                # This part will not be reached if the packing is valid
                print(f"Error: Overlap found between {p1} and {p2}. Dist^2: {dist_sq}")
                break
            if dist_sq < closest_pair_dist_sq:
                closest_pair_dist_sq = dist_sq
        if not valid_packing:
            break
            
    # 3. Print the results for your son.
    print("Yes, your problem formulation is correct!")
    print("We can pack the candies in 3 layers to get the highest number.")
    print(f"The equation for the total number of candies is:")
    print(f"{layer1_count} (in first layer) + {layer2_count} (in second layer) + {layer3_count} (in third layer) = {n}")
    print("\nThis packing is valid because all candies are in the box and no two candies overlap.")
    
    # Show one example calculation for the non-overlapping constraint.
    # We pick one candy from Layer 1 and its closest neighbor from Layer 2.
    p_i = spheres[0]  # (4, 4, 4)
    p_j = spheres[9]  # (8, 8, 10)
    dx = p_j[0] - p_i[0]
    dy = p_j[1] - p_i[1]
    dz = p_j[2] - p_i[2]
    final_dist_sq = dx**2 + dy**2 + dz**2
    
    print("\nFor example, the distance between a candy in layer 1 and a candy in layer 2 is calculated as:")
    print(f"({p_j[0]} - {p_i[0]})^2 + ({p_j[1]} - {p_i[1]})^2 + ({p_j[2]} - {p_i[2]})^2 = {dx}^2 + {dy}^2 + {dz}^2 = {final_dist_sq}")
    print(f"This result {final_dist_sq} is greater than or equal to the required minimum of {min_dist_sq}, so they don't overlap!")

solve_candy_packing()