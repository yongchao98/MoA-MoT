def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes made by a
    majority voting algorithm in the described experts problem.

    Args:
      n (int): The total number of experts.
      c (int): The number of mistakes an expert can make before being removed.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 1 or c < 1:
        print("n and c must be integers greater than or equal to 1.")
        return

    # Bound for M_W: mistakes when the true expert is wrong.
    # M_W <= c - 1
    bound_mw = c - 1

    # Bound for M_R: mistakes when the true expert is right.
    # M_R <= (n - 1) * c / 2
    bound_mr = (n - 1) * c / 2

    # Total bound M = M_W + M_R
    total_bound = bound_mw + bound_mr

    print("Derivation of the upper bound M:")
    print("M <= M_W + M_R")
    print("-" * 20)
    print(f"The number of algorithm mistakes when the true expert is wrong (M_W) is at most c - 1.")
    print(f"M_W <= {c} - 1 = {bound_mw}")
    print("-" * 20)
    print(f"The number of algorithm mistakes when the true expert is right (M_R) is at most (n - 1) * c / 2.")
    print(f"M_R <= ({n} - 1) * {c} / 2 = {bound_mr}")
    print("-" * 20)
    print(f"The total upper bound on mistakes is:")
    print(f"M <= {bound_mw} + {bound_mr} = {total_bound}")

# --- Example Usage ---
# You can change these values to see the bound for different scenarios.
# n: total number of experts
# c: mistake threshold for removing an expert
example_n = 21
example_c = 10

calculate_mistake_bound(example_n, example_c)