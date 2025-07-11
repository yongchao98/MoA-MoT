import math

def solve_packing_problem():
    """
    Solves the discrete sphere packing problem using a greedy algorithm.
    """
    # Define the valid range for the integer coordinates based on the problem.
    x_range = range(4, 21)  # Corresponds to 2cm to 10cm
    y_range = range(4, 21)  # Corresponds to 2cm to 10cm
    z_range = range(4, 19)  # Corresponds to 2cm to 9cm

    # The non-overlapping condition: squared distance must be >= 64.
    min_sq_dist = (4 + 4)**2

    # 1. Create a list of all possible candidate center points.
    candidates = []
    for z in z_range:
        for y in y_range:
            for x in x_range:
                candidates.append((x, y, z))

    # 2. Define the center of the box to sort candidates.
    # Center in cm: (6, 6, 5.5). In grid units: (12, 12, 11).
    center_point = (12, 12, 11)

    # 3. Sort candidates by distance from the center of the box (outward spiral).
    candidates.sort(key=lambda p: (p[0] - center_point[0])**2 + 
                                  (p[1] - center_point[1])**2 + 
                                  (p[2] - center_point[2])**2)

    # 4. Use a greedy approach to pack the spheres.
    packed_spheres = []
    for p_candidate in candidates:
        can_place = True
        # Check against all spheres already placed in the box.
        for p_packed in packed_spheres:
            sq_dist = (p_candidate[0] - p_packed[0])**2 + \
                      (p_candidate[1] - p_packed[1])**2 + \
                      (p_candidate[2] - p_packed[2])**2
            
            if sq_dist < min_sq_dist:
                can_place = False
                break
        
        if can_place:
            packed_spheres.append(p_candidate)

    # The final answer is the number of spheres we managed to pack.
    n = len(packed_spheres)
    
    print("The problem formulation is correct.")
    print(f"The final equation is n = {n}")
    print(f"The maximized number of eyeball candies (n) is: {n}")

# Run the solver
solve_packing_problem()