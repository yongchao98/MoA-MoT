def solve():
    """
    This function calculates and prints the upper bound for d for several values of n
    to illustrate the derivation.
    """
    print("The condition is that for any n, the segments created by {a_1, ..., a_n}")
    print("on [0, d] are all of length at most 1/n.")
    print("This implies d <= (n+1)/n = 1 + 1/n.")
    print("This must hold for all n.\n")

    # Print the bound for some example values of n
    example_n_values = [1, 2, 3, 5, 10, 20, 50, 100]
    for n in example_n_values:
        bound = (n + 1) / n
        print(f"For n = {n:3}, d must be <= ({n}+1)/{n} = {1 + 1/n:.4f}")

    print("\nAs n approaches infinity, the term 1/n approaches 0.")
    print("The limit of the upper bound is 1.")
    print("Therefore, the largest possible value for d is 1.")

solve()