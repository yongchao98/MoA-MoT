import math

def solve_problem():
    """
    Solves the given mathematical problem by analyzing its components and
    performing the necessary calculations.
    """

    # The problem asks to calculate the value of (Sum) * (Integral).

    # Part 1: The Summation
    # The sum is given by S = sum_{i=1 to 10} sum_{j=1 to 10} l(p_(21367+i), p_(14567+j)).
    # The manifold M(n,p) is the Stiefel manifold of orthonormal p-frames in R^n.
    # Its injectivity radius, l(n,p), is a standard result from differential geometry.
    # For n > p, the injectivity radius is pi.
    # In the sum, n = p_(21367+i) and p = p_(14567+j). Since 21367+i > 14567+j for all i,j in [1,10],
    # the condition n > p is always met. The condition n,p >= 5 is also met as the prime indices are large.
    # Therefore, l(n,p) = pi for every term in the sum.
    # The sum simplifies to S = 10 * 10 * pi = 100 * pi.
    sum_factor = 100
    sum_value = sum_factor * math.pi

    # Part 2: The Integral
    # The integral's value depends on dimensions d1 and d2. Let's calculate them.
    # This requires finding prime numbers for large indices.

    def sieve(limit):
        """Generates a list of prime numbers up to a given limit."""
        primes = [True] * (limit + 1)
        if limit >= 0:
            primes[0] = False
        if limit >= 1:
            primes[1] = False
        for i in range(2, int(math.sqrt(limit)) + 1):
            if primes[i]:
                for multiple in range(i * i, limit + 1, i):
                    primes[multiple] = False
        prime_numbers = [i for i, is_prime in enumerate(primes) if is_prime]
        return prime_numbers

    # We need primes up to index 10231. A safe upper bound for the sieve is 280,000.
    # p_n ~ n(log n + log log n) => p_10231 ~ 10231 * (ln(10231)+ln(ln(10231))) ~ 107340
    # Let's use a larger sieve limit just in case.
    prime_list = sieve(300000)

    def get_prime(k):
        """Returns the k-th prime number (1-indexed)."""
        return prime_list[k - 1]

    def dim_M(n, p):
        """Calculates the dimension of the Stiefel manifold M(n,p)."""
        return n * p - p * (p + 1) // 2

    # Calculate d1 and d2
    n1_idx, p1_idx = 8231, 781
    n1, p1 = get_prime(n1_idx), get_prime(p1_idx)
    d1 = dim_M(n1, p1)

    n2_idx, p2_idx = 10231, 2321
    n2, p2 = get_prime(n2_idx), get_prime(p2_idx)
    d2 = dim_M(n2, p2)

    # The integral can be split into two parts by separating the numerator:
    # I = Integral[ Part1 + Part2 ]dx
    # Part2 simplifies to x*exp(-x). The integral of Part2 from 0 to infinity is Gamma(2) = 1.
    # Part1 involves the term (x^(2d1) - x^(2d2)) / (1+x^(2d1))(1+x^(2d2)), which can be rewritten
    # as 1/(1+x^(2d2)) - 1/(1+x^(2d1)).
    # For very large d, the function f(x,d) = 1/(1+x^(2d)) acts like a step function, being ~1 for x<1
    # and ~0 for x>1. Since d1 and d2 are enormous, the difference between these two step-like
    # functions is effectively zero everywhere.
    # Thus, the integral of Part1 is 0.
    # The total value of the integral is I = 0 + 1 = 1.
    integral_value = 1.0

    # Part 3: Final Calculation
    final_result = sum_value * integral_value
    
    print("This problem simplifies as follows:")
    print(f"1. The sum term evaluates to {sum_factor} * pi, as the injectivity radius l(n,p) is pi.")
    print(f"2. The integral term evaluates to {integral_value}.")
    print("   - The dimensions in the integrand are d1 = {d1} and d2 = {d2}.".format(d1=d1, d2=d2))
    print("   - Due to these large exponents, one part of the integral becomes 0.")
    print("   - The other part, Integral[x*exp(-x)], is exactly 1.")
    
    print("\nThe final calculation is:")
    print(f"({sum_factor} * {math.pi}) * {integral_value} = {final_result}")

solve_problem()