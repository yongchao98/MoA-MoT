def solve_experts_problem_bound():
    """
    Calculates and prints an upper bound on the number of mistakes made by a
    majority voting algorithm in a variant of the experts problem.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert can make before being removed.
    """
    # Parameters for the problem, you can change these values
    n = 10
    c = 5

    # An expert is removed after making c mistakes.
    # A true expert makes strictly fewer than c mistakes, i.e., at most c-1.

    # Bound for mistakes when the true expert is wrong.
    # This is bounded by the max number of mistakes the true expert can make.
    mw_bound = c - 1

    # Bound for mistakes when the true expert is correct but outvoted.
    # This is bounded by the total mistake "budget" of all n-1 false experts.
    mc_bound = c * (n - 1)

    # The total upper bound is the sum of the two bounds.
    total_bound = mw_bound + mc_bound

    print(f"Given n={n} experts and a mistake threshold c={c}:")
    print("An upper bound on the number of mistakes M is the sum of:")
    print(" 1. The maximum mistakes when the true expert is wrong (M_W).")
    print(" 2. The maximum mistakes when the true expert is correct but outvoted (M_C).")
    print("")
    # Print the final equation with each number explicitly shown
    print("M <= M_W + M_C")
    print(f"M <= ({c} - 1) + {c} * ({n} - 1)")
    print(f"M <= {mw_bound} + {mc_bound}")
    print(f"M <= {total_bound}")

solve_experts_problem_bound()
