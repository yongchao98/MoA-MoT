def solve_divisor_problem():
    """
    This function calculates the size of the largest union of 20 antichains
    in the divisor poset of N = 823564528378596.

    The problem is interpreted as operating on the set of divisors of N.
    The key steps are:
    1.  By the dual of Dilworth's theorem, a union of 20 antichains is a set
        where the longest chain has a length of at most 20.
    2.  The prime factorization of N is 2^2 * (product of 15 other primes).
    3.  The length of the longest chain in the divisor poset of N is Omega(N) + 1,
        where Omega(N) is the sum of the exponents of its prime factors.
        Omega(N) = 2 + 15 = 17. So, the longest chain has length 18.
    4.  Since 18 <= 20, the entire set of divisors of N satisfies the condition.
    5.  The size of this set is tau(N), the number of divisors.
        tau(N) = (2+1) * (1+1)^15 = 3 * 2^15.
    """

    # Exponents from the prime factorization of N = 823564528378596
    # N = 2^2 * 3^1 * 7^1 * ... * 59^1 (15 primes other than 2)
    exponent_of_2 = 2
    num_other_primes = 15
    exponent_of_other_primes = 1

    # Calculate tau(N) = (exponent_of_2 + 1) * (exponent_of_other_primes + 1)^num_other_primes
    factor1 = exponent_of_2 + 1
    base = exponent_of_other_primes + 1
    exponent = num_other_primes

    result = factor1 * (base ** exponent)

    # Output the final equation with each number, as requested.
    print(f"{factor1} * {base}**{exponent} = {result}")

solve_divisor_problem()