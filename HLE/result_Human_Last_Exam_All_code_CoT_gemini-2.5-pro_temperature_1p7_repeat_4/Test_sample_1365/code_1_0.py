def solve_experts_problem(n, c):
    """
    Calculates and prints the upper bound on the number of mistakes for a variant of the experts problem.

    Args:
        n (int): The total number of experts.
        c (int): The mistake threshold for removing an expert.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n <= 1 or c <= 0:
        print("Error: n must be an integer greater than 1, and c must be a positive integer.")
        return

    # Calculate the maximum total number of mistakes that all experts can make combined.
    # W_max = (n-1)*c (for non-true experts) + (c-1) (for the true expert)
    max_total_expert_mistakes = n * c - 1

    # The number of algorithm mistakes (M) must satisfy M < 2 * W_max.
    # Since M is an integer, the upper bound is 2 * W_max - 1.
    bound = 2 * max_total_expert_mistakes - 1

    # The number of mistakes cannot be negative.
    if bound < 0:
        bound = 0

    print(f"An upper bound on the number of mistakes (M) for n={n} experts and a mistake threshold of c={c} can be calculated as follows:")
    print()
    print("M < 2 * (Total Expert Mistake Budget)")
    print("M < 2 * ((n-1)*c + (c-1))")
    print("M < 2 * (n*c - 1)")
    print()
    print("Solving the final equation for the integer bound:")
    print(f"M <= 2 * (n * c - 1) - 1")
    print(f"M <= 2 * ({n} * {c} - 1) - 1")
    print(f"M <= 2 * ({n * c} - 1) - 1")
    print(f"M <= 2 * {max_total_expert_mistakes} - 1")
    print(f"M <= {2 * max_total_expert_mistakes} - 1")
    print(f"M <= {bound}")

# You can change these values to explore the bound with different parameters.
n_experts = 10
c_mistakes = 5

solve_experts_problem(n_experts, c_mistakes)
<<<97>>>