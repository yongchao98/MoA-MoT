def solve():
    """
    Calculates the number of internal adjunctions in the given 2-category
    from F_11^3 to itself.

    This problem reduces to calculating the number of automorphisms of the
    vector space F_11^3, which is the order of the general linear group GL(3, F_11).
    The formula for the order of GL(n, q) is:
    |GL(n, q)| = (q^n - 1)(q^n - q)(q^n - q^2)...(q^n - q^(n-1))
    """
    q = 11
    n = 3

    # Calculate the terms of the formula
    # |GL(3, 11)| = (11^3 - 11^0) * (11^3 - 11^1) * (11^3 - 11^2)
    q_n = q**n
    
    term1 = q_n - (q**0)
    term2 = q_n - (q**1)
    term3 = q_n - (q**2)
    
    # Calculate the final result
    result = term1 * term2 * term3
    
    # Print the result in the specified equation format
    print(f"{term1} * {term2} * {term3} = {result}")

solve()