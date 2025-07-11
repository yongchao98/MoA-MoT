def solve_for_max_n():
    """
    Calculates the maximum value of n for a given k based on the formula
    n = k^2 - k + 1.
    """
    # The user can change this value to see the result for a different k.
    # We assume k >= 2 as is standard for such problems.
    k = 5

    print(f"The problem asks for the maximum value of n in terms of k for a k-uniform intersecting family with full differences of size k-1.")
    print(f"The result from extremal set theory is that the maximum value is n = k^2 - k + 1.")
    print(f"Here is the calculation for the example case k = {k}:\n")

    # Perform the calculation
    k_squared = k**2
    n = k_squared - k + 1

    # Output the equation with each number, as requested.
    print("The final equation is:")
    print(f"n = {k}^2 - {k} + 1")
    print(f"n = {k_squared} - {k} + 1")
    print(f"n = {n}")

solve_for_max_n()