import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary {prime: exponent}.
    """
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def solve():
    """
    Calculates the largest number of mutually independent events for rolling
    100 6-sided dice.
    """
    base = 6
    power = 100

    # The size of the sample space is base^power.
    # The maximum number of mutually independent events is the sum of the exponents
    # in the prime factorization of the sample space size.
    # N = (p1^a1 * p2^a2 * ...)^power = p1^(a1*power) * p2^(a2*power) * ...
    # The result is sum(ai * power).

    # First, find the prime factorization of the base.
    base_factors = get_prime_factorization(base)

    # The exponents in the prime factorization of base^power are the exponents
    # of the base's factors, each multiplied by the power.
    # The largest number of events 'm' is the sum of these new exponents.
    m = 0
    
    # Store calculation details for the explanation
    parts = []
    
    for prime, exponent in base_factors.items():
        total_exponent = exponent * power
        m += total_exponent
        parts.append(str(total_exponent))
        
    print(f"The size of the sample space is {base}^{power}.")
    print(f"The prime factorization of the base ({base}) is {base_factors}.")
    print(f"The prime factorization of the sample space size ({base}^{power}) has exponents given by the exponents of the base's factors multiplied by {power}.")
    print(f"The largest possible value of m is the sum of these exponents:")
    print(f"m = {' + '.join(parts)} = {m}")


solve()