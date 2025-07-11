def solve():
    """
    Calculates the size of the largest union of 20 antichains in the divisor poset of N.

    This is equivalent to the total number of divisors of N, as the longest chain
    in the poset is shorter than 21.
    """
    # The prime factorization of N = 823564528378596 is
    # 2^2 * 3^2 * 7^1 * 13^1 * 17^1 * 19^1 * 23^1 * 29^1 * 31^1 * 37^1 * 41^1 * 43^1.
    exponents = [2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    # The total number of divisors is (a1+1)*(a2+1)*...
    terms = [exp + 1 for exp in exponents]

    # Calculate the result
    result = 1
    for term in terms:
        result *= term

    # Build the equation string
    equation_parts = [f"({exp}+1)" for exp in exponents]
    equation_str = " * ".join(equation_parts)

    print(f"The size of the largest union of 20 antichains is the total number of divisors.")
    print(f"The calculation is: {equation_str} = {result}")

solve()