def solve_compactness_problem():
    """
    Calculates the n-compactness number for X = [0,1]^3.

    The value [X] is the subbase compactness number, s(X).
    For a product of k compact Hausdorff spaces, the formula is:
    s(X_1 x ... x X_k) = sum(s(X_i)) - k + 1

    For the interval I = [0,1], s(I) = 2.
    For X = [0,1]^3, k=3 and each s(X_i) = 2.
    """

    # Number of spaces in the product
    k = 3

    # Subbase compactness number for the interval [0,1]
    s_I = 2

    # The list of s(X_i) values
    s_values = [s_I] * k

    # Sum of s(X_i)
    sum_s_values = sum(s_values)

    # Apply the formula: s(X) = sum(s(X_i)) - k + 1
    result = sum_s_values - k + 1

    # Print the final equation with all numbers
    equation_parts = [str(s) for s in s_values]
    print(f"The n-compactness number [X] for X = [0,1]^3 is calculated using the product formula.")
    print(f"First, we note that for the interval I = [0,1], the value [I] = {s_I}.")
    print(f"For the product space X = [0,1]^3, the formula gives:")
    print(f"[X] = {' + '.join(equation_parts)} - {k} + 1 = {result}")

solve_compactness_problem()