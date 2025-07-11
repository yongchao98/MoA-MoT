import math
import sympy

def prime_sieve(n):
    """Generates primes up to n using a sieve."""
    primes = [True] * (n + 1)
    primes[0] = primes[1] = False
    for i in range(2, int(math.sqrt(n)) + 1):
        if primes[i]:
            for multiple in range(i*i, n + 1, i):
                primes[multiple] = False
    prime_numbers = [i for i, is_p in enumerate(primes) if is_p]
    return prime_numbers

def get_ith_prime(i):
    """
    Gets the i-th prime number.
    Uses sympy.prime for simplicity and accuracy.
    A manual implementation would require a large sieve or estimations.
    """
    return sympy.prime(i)

def get_dim(n, p):
    """Calculates the dimension of the manifold M(n,p)."""
    return n * p - p * (p + 1) // 2

def solve():
    """
    Solves the given mathematical problem.
    """
    # Term 1: The Summation
    # The term l(n,p) is the injectivity radius of the Stiefel manifold M(n,p).
    # This is a standard result in differential geometry, and its value is pi.
    # It is independent of n and p (for n > p >= 1).
    # The condition n, p >= 5 is satisfied by the primes used.
    l_val = math.pi
    sum_term = 0
    # The problem asks for the sum over i=1..10 and j=1..10.
    # Since l(n,p) is constant (pi), the sum is 10 * 10 * pi.
    sum_term = 10 * 10 * l_val

    print(f"Step 1: Calculate the summation term.")
    print(f"The term l(n,p) represents the injectivity radius of the Stiefel manifold M(n,p), which is pi.")
    print(f"The summation is sum_{{i=1}}^{{10}} sum_{{j=1}}^{{10}} l(n,p) = 10 * 10 * pi = 100 * pi.")
    print(f"Value of the summation term = {sum_term}\n")

    # Term 2: The Integral
    # The integral can be split into two parts.
    # Part A: integral of x * exp(-x) from 0 to infinity
    integral_part_A = 1.0  # This is Gamma(2) = 1! = 1

    # Part B: The complex fraction involving dimensions d1 and d2.
    # Let's calculate d1 and d2.
    p_8231 = get_ith_prime(8231)
    p_781 = get_ith_prime(781)
    n1 = p_8231
    p1 = p_781
    d1 = get_dim(n1, p1)

    p_10231 = get_ith_prime(10231)
    p_2321 = get_ith_prime(2321)
    n2 = p_10231
    p2 = p_2321
    d2 = get_dim(n2, p2)

    print(f"Step 2: Calculate the integral term.")
    print("The integral splits into two parts: I = I1 + I2.")
    print(f"I2 = integral from 0 to inf of x*exp(-x) dx = 1.\n")
    print("I1 is the integral of the first fraction.")
    print("Let's calculate the dimensions d1 and d2 involved in I1.")
    print(f"n1 = p_8231 = {n1}")
    print(f"p1 = p_781 = {p1}")
    print(f"d1 = dim(M(n1,p1)) = n1*p1 - p1*(p1+1)/2 = {d1}")
    print(f"n2 = p_10231 = {n2}")
    print(f"p2 = p_2321 = {p2}")
    print(f"d2 = dim(M(n2,p2)) = n2*p2 - p2*(p2+1)/2 = {d2}\n")

    print("The dimensions d1 and d2 are extremely large.")
    print("The integrand of I1 contains the term (x^(2d1) - x^(2d2)) / ((1+x^(2d1))*(1+x^(2d2))).")
    print("As d1 and d2 approach infinity, this term approaches 0 for all x > 0.")
    print("By the Dominated Convergence Theorem, the integral I1 evaluates to 0.")
    integral_part_B = 0.0

    integral_term = integral_part_B + integral_part_A
    print(f"So, the total value of the integral is I1 + I2 = 0 + 1 = {integral_term}\n")

    # Final Result
    final_result = sum_term * integral_term
    print(f"Step 3: Calculate the final result.")
    print("The final result is the product of the summation term and the integral term.")
    print(f"Final Result = (100 * pi) * 1")
    print(f"Final equation: ({100} * {math.pi}) * ({integral_part_B} + {integral_part_A}) = {final_result}")

solve()
<<<314.1592653589793>>>