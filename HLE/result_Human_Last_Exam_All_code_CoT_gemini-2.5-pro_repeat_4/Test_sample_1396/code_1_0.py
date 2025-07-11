def solve_cake_cutting_bound():
    """
    Calculates the upper bound on query complexity for envy-free cake cutting
    for 4 agents based on the CGMM22 protocol.
    """
    # Number of agents
    n = 4

    # The number of 'cut' queries is given by the formula C(n) = n^2 - 1.
    num_cuts = n**2 - 1

    # The number of 'eval' queries is given by the formula E(n) = n * (n^2 - 1) / 3.
    # The result is guaranteed to be an integer for any integer n.
    num_evals = int(n * (n**2 - 1) / 3)

    # The total upper bound O is the sum of the number of cut and eval queries.
    total_queries = num_cuts + num_evals

    # Print the final equation with all the numbers, as requested.
    print(f"The upper bound O for n={n} agents is based on the state-of-the-art O(n^2) protocol.")
    print("The calculation is as follows:")
    print(f"Number of 'cut' queries = {n}^2 - 1 = {num_cuts}")
    print(f"Number of 'eval' queries = {n} * ({n}^2 - 1) / 3 = {num_evals}")
    print(f"Total Upper Bound O = {num_cuts} + {num_evals} = {total_queries}")

solve_cake_cutting_bound()