import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B with n sides of lengths a_i.

    Args:
        n (int): The number of sides of the polygon B.
        a (list of float): The lengths of the edges of B, [a_1, a_2, ..., a_n].
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return
    if len(a) != n:
        print(f"The number of edge lengths ({len(a)}) does not match n ({n}).")
        return

    # phi is the angle between consecutive normals
    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    cos_phi_half = math.cos(phi / 2)
    tan_phi_quarter = math.tan(phi / 4)

    # Use a 0-based list for 'a', but calculations are 1-based as in the formula.
    # We create a temporary 1-based list a_1based to match the math notation.
    a_1based = [0] + a

    distances = []
    print("Calculating the distance component for each vertex:")

    for i in range(1, n + 1):
        # Handle wraparound for indices (n+1 becomes 1)
        i_plus_1 = i + 1 if i < n else 1

        ai = a_1based[i]
        ai_plus_1 = a_1based[i_plus_1]
        
        # Calculate b_i, the length of the diagonal skipping vertex i.
        # The formula uses the index convention where b_i depends on a_i and a_{i+1}.
        # My vertex indexing would make b_i depend on a_i and a_{i+1}. Let's assume the problem
        # intended a_i and a_{i+1} to be the adjacent edges.
        b_i_sq = ai**2 + ai_plus_1**2 + 2 * ai * ai_plus_1 * cos_phi
        b_i = math.sqrt(b_i_sq)

        # Calculate the Hausdorff distance component D_i
        numerator = ai * ai_plus_1 * cos_phi_half * tan_phi_quarter
        d_i = numerator / b_i
        distances.append(d_i)

        # Print the calculation for the current D_i
        print(f"D_{i} = ({ai} * {ai_plus_1} * cos({phi/2:.4f}) * tan({phi/4:.4f})) / sqrt({ai}^2 + {ai_plus_1}^2 + 2*{ai}*{ai_plus_1}*cos({phi:.4f}))")
        print(f"    = ({ai} * {ai_plus_1} * {cos_phi_half:.4f} * {tan_phi_quarter:.4f}) / {b_i:.4f}")
        print(f"    = {d_i:.6f}")

    max_dist = max(distances)
    print("\nThe largest possible Hausdorff distance is the maximum of these values.")
    print(f"Max Hausdorff Distance = {max_dist:.6f}")
    return max_dist

if __name__ == '__main__':
    # Example: A rectangle with sides 3 and 4
    n_sides = 4
    edge_lengths = [3.0, 4.0, 3.0, 4.0]
    print(f"Running example for n={n_sides} with edge lengths a={edge_lengths}\n")
    calculate_max_hausdorff_distance(n_sides, edge_lengths)
    
    print("\n" + "="*50 + "\n")

    # Example: A regular triangle (n=3) with side length 5
    n_sides_triangle = 3
    edge_lengths_triangle = [5.0, 5.0, 5.0]
    print(f"Running example for n={n_sides_triangle} with edge lengths a={edge_lengths_triangle}\n")
    calculate_max_hausdorff_distance(n_sides_triangle, edge_lengths_triangle)