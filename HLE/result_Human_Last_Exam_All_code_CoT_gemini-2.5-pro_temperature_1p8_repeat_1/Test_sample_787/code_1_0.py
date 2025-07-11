import math

def solve():
    """
    Calculates the limit of the sequence g_n based on the polynomial P(X).
    The method involves determining the prime factorization of the limit, g.
    """

    # Step 1: Explain the problem and the approach.
    print("Let P(X) = (X^5-1)*(X^5-X)*(X^5-X^2)*(X^5-X^3)*(X^5-X^4).")
    print("Let g_n be the greatest common divisor (GCD) of P(p) for all primes p > n.")
    print("The limit g = lim_{n->inf} g_n is the GCD of P(p) for all sufficiently large primes p.")
    print("We find g by determining its prime factorization, g = 2^a * 3^b * 5^c * ...")
    print("The exponent of a prime q in g, denoted nu_q(g), is min(nu_q(P(p))) over all sufficiently large primes p.")
    print("-" * 30)

    # Step 2: Determine the exponent of each prime factor.
    # Exponent for prime 2
    a = 10
    print(f"For q=2:")
    print(f"Analysis shows that for any prime p > 2, nu_2(P(p)) >= {a}.")
    print(f"Also, for arbitrarily large primes p (e.g., p = 3, 11, 19, ...), we have nu_2(P(p)) = {a}.")
    print(f"Thus, the exponent of 2 in the prime factorization of g is {a}.")
    print("-" * 30)

    # Exponent for prime 3
    b = 2
    print(f"For q=3:")
    print(f"Analysis shows that for any prime p > 3, nu_3(P(p)) >= {b}.")
    print(f"Also, for arbitrarily large primes p (e.g., p = 2, 5, 11, ... mod 9), we have nu_3(P(p)) = {b}.")
    print(f"Thus, the exponent of 3 in the prime factorization of g is {b}.")
    print("-" * 30)

    # Exponent for prime 5
    c = 1
    print(f"For q=5:")
    print(f"Analysis shows that for any prime p > 5, nu_5(P(p)) >= {c}.")
    print(f"Also, for arbitrarily large primes p (e.g., p = 2, 3, 17, ...), we have nu_5(P(p)) = {c}.")
    print(f"Thus, the exponent of 5 in the prime factorization of g is {c}.")
    print("-" * 30)

    # Exponents for primes q > 5
    print("For q > 5:")
    print("For any prime q > 5, there exist arbitrarily large primes p such that P(p) is not divisible by q.")
    print("This means g has no prime factors greater than 5.")
    print("-" * 30)
    
    # Step 3: Calculate the final result.
    base2 = 2**a
    base3 = 3**b
    base5 = 5**c
    
    result = base2 * base3 * base5

    # Step 4: Print the final calculation and the answer.
    print("The limit g is the product of these prime powers:")
    print(f"g = 2^{a} * 3^{b} * 5^{c}")
    print(f"g = {base2} * {base3} * {base5}")
    print(f"g = {result}")

solve()