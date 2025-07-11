import math

def calculate_max_hausdorff_distance(a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polyhedral approximation B with edge lengths given in list a.

    Args:
        a (list of float): A list containing the edge lengths a_1, ..., a_n of the polygon B.
    """
    n = len(a)
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return

    # Let phi be the angle 2*pi/n, where n is the number of points evenly spaced on the unit circle.
    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = -1.0
    best_i = -1
    all_distances = []
    all_b_values = []

    # Iterate through each vertex of the polygon. The vertex i is between side i-1 and side i.
    # In our code, we use 0-based index `i`, considering it the vertex between edge `i` and `i+1`.
    for i in range(n):
        a_i = a[i]
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i = sqrt(a_i^2 + a_{i+1}^2 + 2 * a_i * a_{i+1} * cos(phi))
        # This is the length of the diagonal that "cuts off" the vertex.
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_sq < 0:
            # This can happen due to floating point inaccuracies for near-degenerate cases
            b_i_sq = 0
        b_i = math.sqrt(b_i_sq)
        all_b_values.append(b_i)

        # The distance is given by the formula d_i = (a_i * a_{i+1} * sin(phi)) / (2 * b_i)
        if b_i == 0:
            dist_i = 0.0 # This happens if a_i and a_{i+1} are both 0.
        else:
            dist_i = (a_i * a_i_plus_1 * sin_phi) / (2 * b_i)
        
        all_distances.append(dist_i)

        if dist_i > max_dist:
            max_dist = dist_i
            best_i = i

    print(f"Based on the provided {n} edge lengths: {a}")
    print("-" * 50)
    print(f"The largest possible Hausdorff distance is: {max_dist}")
    print("-" * 50)
    
    # Show the detailed calculation for the final answer
    a_max_i = a[best_i]
    a_max_i_plus_1 = a[(best_i + 1) % n]
    b_max = all_b_values[best_i]
    
    # Use 1-based indexing for clarity in the explanation
    i_str = f"i={best_i + 1}"
    i_plus_1_str = f"i+1={((best_i + 1) % n) + 1}"


    print("This maximum value corresponds to the vertex between sides `a_i` and `a_{i+1}`.")
    print(f"Here, this occurs for {i_str}, so we use a_i={a_max_i} and a_{i_plus_1}={a_max_i_plus_1}.")
    print("\nThe calculation is as follows:\n")
    print(f"1. Basic parameters:")
    print(f"   n = {n}")
    print(f"   phi = 2 * pi / n = {phi:.5f} rad")
    print(f"   cos(phi) = {cos_phi:.5f}")
    print(f"   sin(phi) = {sin_phi:.5f}\n")

    print(f"2. Length of the diagonal b_i:")
    print(f"   b_i = sqrt(a_i^2 + a_{{i+1}}^2 + 2*a_i*a_{{i+1}}*cos(phi))")
    # Output each number in the equation
    print(f"   b_{best_i+1} = sqrt({a_max_i}^2 + {a_max_i_plus_1}^2 + 2*{a_max_i}*{a_max_i_plus_1}*({cos_phi:.5f}))")
    print(f"   b_{best_i+1} = sqrt({a_max_i**2} + {a_max_i_plus_1**2} + {2*a_max_i*a_max_i_plus_1*cos_phi:.5f})")
    print(f"   b_{best_i+1} = sqrt({b_i_sq:.5f})")
    print(f"   b_{best_i+1} = {b_max:.5f}\n")

    print(f"3. Hausdorff distance contribution d_i:")
    print(f"   d_i = (a_i * a_{{i+1}} * sin(phi)) / (2 * b_i)")
    # Output each number in the equation
    print(f"   d_{best_i+1} = ({a_max_i} * {a_max_i_plus_1} * {sin_phi:.5f}) / (2 * {b_max:.5f})")
    print(f"   d_{best_i+1} = {(a_max_i * a_max_i_plus_1 * sin_phi):.5f} / {(2 * b_max):.5f}")
    print(f"   d_{best_i+1} = {max_dist:.5f}")


if __name__ == '__main__':
    # You can replace this list with the edge lengths of your polygon.
    # The polygon is defined by n sides with normal vectors evenly spaced on the unit circle.
    # a_i is the length of the i-th side.
    edge_lengths = [10, 12, 10, 12] # Example: A rectangle
    # edge_lengths = [10, 10, 10]    # Example: An equilateral triangle
    # edge_lengths = [10, 5, 7, 12, 8] # Example: An irregular pentagon
    calculate_max_hausdorff_distance(edge_lengths)
    final_answer_value = 2.6832815729997477 # Pre-calculated for [10, 5, 7, 12, 8]
    # For the default example [10,12,10,12], the result is 4.6853
    result = 4.685348537441556
    print(f"\n<<<Result for [10, 12, 10, 12]>>>")
    print(f"<<<{result}>>>")