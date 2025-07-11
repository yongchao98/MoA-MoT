import math

def get_prime_factorization(n):
    """
    Returns a dictionary where keys are prime factors and values are their exponents.
    For example, get_prime_factorization(12) returns {2: 2, 3: 1}.
    """
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while (temp_n % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def solve_max_independent_events():
    """
    Calculates the maximum number of mutually independent events for rolling
    100 six-sided dice.
    """
    num_dice = 100
    sides_per_die = 6

    # The size of the sample space is N = sides_per_die ^ num_dice.
    # We need the sum of the exponents in the prime factorization of N.
    # N = (p1^a1 * p2^a2 * ...)^num_dice = p1^(a1*num_dice) * p2^(a2*num_dice) * ...
    # The sum of exponents for N is (a1 + a2 + ...) * num_dice.

    # First, find the prime factorization of the base (sides_per_die).
    base_factors = get_prime_factorization(sides_per_die)

    # Calculate the exponents for N = 6^100.
    exponents_in_N = []
    for prime, exponent_in_base in base_factors.items():
        exponent_in_N = exponent_in_base * num_dice
        exponents_in_N.append(exponent_in_N)

    # The maximum number of events 'm' is the sum of these exponents.
    max_m = sum(exponents_in_N)

    # Print the explanation and the final equation.
    print("The size of the sample space is 6^100.")
    print("The prime factorization of 6^100 is (2 * 3)^100 = 2^100 * 3^100.")
    print("The maximum number of mutually independent events is the sum of the exponents.")
    
    # Construct the equation string: "100 + 100 = 200"
    equation_str = " + ".join(map(str, exponents_in_N))
    print(f"The final calculation is: {equation_str} = {max_m}")

solve_max_independent_events()