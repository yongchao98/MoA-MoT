import math

def solve_hausdorff_distance(n, edge_lengths):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B with given edge lengths.

    Args:
        n (int): The number of sides/vertices of the polygon B.
        edge_lengths (list[float]): A list of the lengths of the edges of B.
    """
    if n != len(edge_lengths):
        print("Error: The number of sides n must be equal to the number of edge lengths provided.")
        return
    if n < 3:
        print("Error: A polygon must have at least 3 sides.")
        return

    # Let the given edge lengths be a_1, ..., a_n. In the code, we use 0-based indexing a_0, ..., a_{n-1}.
    a = edge_lengths
    phi = 2 * math.pi / n

    max_dist = -1.0
    best_i = -1

    for i in range(n):
        a_i = a[i]
        # The vertex is between edge i and edge i+1.
        # Let's adjust indexing to match the formula convention a_i, a_{i+1}
        # Let our loop index i correspond to the vertex between edge i and edge i+1
        a_i_val = a[i]
        a_i_plus_1_val = a[(i + 1) % n]
        
        # Avoid issues with very small or zero edge lengths which can make b_i zero.
        if a_i_val == 0 or a_i_plus_1_val == 0:
            dist = 0
        else:
            # Calculate b_i, the length of the diagonal connecting the endpoints
            # of the edges a_i and a_{i+1} that are not at the common vertex.
            # Using formula b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))
            b_i_sq = a_i_val**2 + a_i_plus_1_val**2 + 2 * a_i_val * a_i_plus_1_val * math.cos(phi)
            if b_i_sq < 0: b_i_sq = 0 # handle potential floating point inaccuracies for small negative values
            b_i = math.sqrt(b_i_sq)
    
            if b_i == 0:
              dist = 0.0
            else:
              # This formula gives the largest possible Hausdorff distance local to a vertex
              dist = (a_i_val * a_i_plus_1_val * math.tan(phi / 2)) / (2 * b_i)

        if dist > max_dist:
            max_dist = dist
            best_i = i

    print(f"The number of sides n = {n}")
    print(f"The edge lengths are: {edge_lengths}")
    print(f"The external angle phi = 2*pi/n = {phi:.4f} radians")
    
    a_i = a[best_i]
    a_i_plus_1 = a[(best_i + 1) % n]
    b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * math.cos(phi)
    if b_i_sq < 0: b_i_sq = 0
    b_i = math.sqrt(b_i_sq)

    print("\nThe largest possible Hausdorff distance is maximized at the vertex between edges of")
    print(f"length a_{best_i+1} = {a_i} and a_{best_i+2} = {a_i_plus_1}")
    print("\nThe calculation for the maximum distance is:")
    
    print(f"({a_i} * {a_i_plus_1} * tan(({phi:.4f}) / 2)) / (2 * {b_i:.4f})")
    print(f"= ({a_i * a_i_plus_1:.4f} * {math.tan(phi / 2):.4f}) / {2 * b_i:.4f}")
    print(f"= {a_i * a_i_plus_1 * math.tan(phi / 2):.4f} / {2 * b_i:.4f}")
    print(f"= {max_dist:.4f}")

    print("\nFinal Answer:")
    # Using format specifier to avoid trailing zeros
    print(f"{max_dist:.10f}".rstrip('0').rstrip('.'))
    
    # Returning the final answer in the requested format
    return f"<<<{max_dist}>>>"

# Example Usage:
# You can change these values to solve for your specific case.
n_sides = 5
# edge_lengths = [1.0, 1.0, 1.0, 1.0, 1.0] # Regular pentagon
edge_lengths = [3.0, 4.0, 5.0, 2.0, 6.0]   # An irregular pentagon

solve_hausdorff_distance(n_sides, edge_lengths)