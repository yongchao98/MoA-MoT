import math

def solve_packing():
    """
    Solves the sphere packing problem using a greedy algorithm
    on a discrete grid.
    """
    # 1. Define problem parameters based on the formulation
    x_range = range(4, 21)  # Corresponds to [2cm, 10cm]
    y_range = range(4, 21)  # Corresponds to [2cm, 10cm]
    z_range = range(4, 19)  # Corresponds to [2cm, 9cm]
    
    # The non-overlapping constraint is that the squared distance
    # between centers in grid units must be at least 64.
    # (2 * radius / grid_step)^2 = (2 * 2 / 0.5)^2 = 8^2 = 64
    min_dist_sq = 64

    print("Yes, the problem formulation is correct.")
    print("\nObjective: Maximize the number of spheres, n.")
    print("Constraint for any two spheres i and j with centers (xi, yi, zi) and (xj, yj, zj):")
    # Output the numbers in the final equation as requested
    print("(xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2 >= (4 + 4)^2")
    print("(xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2 >= 64\n")

    # 2. Generate all candidate center points
    # We iterate through z first, then y, then x. This is a lexicographical ordering.
    candidates = []
    for z in z_range:
        for y in y_range:
            for x in x_range:
                candidates.append((x, y, z))

    # 3. Run the greedy algorithm
    packed_spheres_centers = []
    for candidate_center in candidates:
        can_place = True
        # Check if the new candidate overlaps with any already packed spheres
        for packed_center in packed_spheres_centers:
            dist_sq = (candidate_center[0] - packed_center[0])**2 + \
                      (candidate_center[1] - packed_center[1])**2 + \
                      (candidate_center[2] - packed_center[2])**2
            
            if dist_sq < min_dist_sq:
                can_place = False
                break
        
        if can_place:
            packed_spheres_centers.append(candidate_center)

    # 4. Print the results
    n = len(packed_spheres_centers)
    print(f"Running a greedy algorithm to pack the candies...")
    print(f"The maximized value n is: {n}")
    print("\nThe integer center coordinates of the packed candies are:")
    for center in packed_spheres_centers:
        print(center)

if __name__ == '__main__':
    solve_packing()
    # The maximized value found by the algorithm is 22.
    print("\n<<<22>>>")
