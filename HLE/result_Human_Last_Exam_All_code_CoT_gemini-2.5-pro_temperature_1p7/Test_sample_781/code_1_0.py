import math

def solve_continuum_problem():
    """
    Calculates the largest number n for the described continuum decomposition problem.
    The problem is a known result in continuum theory. Given a continuum X
    and a set of 'm' distinct points such that no proper subcontinuum contains
    any 'k' of them, the largest number 'n' of continua in an irreducible cover
    X = A_1 U ... U A_n is given by the binomial coefficient C(m-1, k-1).

    In this problem:
    m = 5 (the number of points a, b, c, d, e)
    k = 3 (the size of the subsets of points not contained in any proper subcontinuum)
    """

    m = 5
    k = 3

    # The formula is n = C(m-1, k-1)
    m_prime = m - 1
    k_prime = k - 1
    
    # Calculate the binomial coefficient C(n, k) = n! / (k! * (n-k)!)
    n = math.comb(m_prime, k_prime)

    print(f"The number of special points is m = {m}.")
    print(f"The size of the point subsets defining the property is k = {k}.")
    print(f"The largest number n is given by the formula C(m-1, k-1).")
    print(f"n = C({m}-1, {k}-1) = C({m_prime}, {k_prime})")
    print(f"The calculation is: {m_prime}! / ({k_prime}! * ({m_prime}-{k_prime})!)")
    
    result_str = f"{math.factorial(m_prime)} / ({math.factorial(k_prime)} * {math.factorial(m_prime-k_prime)}) = {n}"
    
    print(f"The final equation is: 4! / (2! * 2!) = {result_str}")
    print(f"The largest number n is: {n}")

solve_continuum_problem()