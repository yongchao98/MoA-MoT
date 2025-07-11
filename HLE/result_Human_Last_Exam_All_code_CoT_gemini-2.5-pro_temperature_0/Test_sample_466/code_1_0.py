def solve():
    """
    Calculates the number of internal adjunctions in the given 2-category
    from F_11^3 to itself. This number is equal to the size of the
    general linear group GL(3, F_11).
    """
    q = 11
    n = 3

    # Calculate the terms in the formula for the size of GL(n, q)
    q_n = q**n
    term1 = q_n - q**0
    term2 = q_n - q**1
    term3 = q_n - q**2

    # Calculate the final result
    result = term1 * term2 * term3

    # Print the equation with all numbers
    print(f"The number of adjunctions is the size of GL({n}, {q}), calculated as:")
    print(f"({q}^{n} - {q}^0) * ({q}^{n} - {q}^1) * ({q}^{n} - {q}^2) = "
          f"({q_n} - {q**0}) * ({q_n} - {q**1}) * ({q_n} - {q**2})")
    print(f"= {term1} * {term2} * {term3}")
    print(f"= {result}")

solve()