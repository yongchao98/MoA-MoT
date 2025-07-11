import math

def get_prime_factorization_exponents(n):
    """
    Finds the exponents of the prime factors of a number n.
    Returns a dictionary of {prime: exponent}.
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

def solve_max_independent_events():
    """
    Calculates the largest possible number of mutually independent events
    for rolling a number of n-sided dice.
    """
    num_dice = 100
    num_sides = 6

    # Get the prime factorization of the number of sides.
    side_factors = get_prime_factorization_exponents(num_sides)

    # The exponents in the prime factorization of |Omega| = num_sides**num_dice
    # are the exponents of the factors of num_sides multiplied by num_dice.
    total_exponents = []
    for prime in sorted(side_factors.keys()):
        exponent = side_factors[prime]
        total_exponent_for_prime = exponent * num_dice
        total_exponents.append(total_exponent_for_prime)

    # The maximum number of events is the sum of these total exponents.
    max_m = sum(total_exponents)

    # Format the output to show the final equation as requested.
    equation_parts = [str(e) for e in total_exponents]
    print(f"The size of the sample space is {num_sides}^{num_dice}.")
    print(f"The prime factorization of {num_sides} is {' * '.join([f'{p}^{e}' for p, e in side_factors.items()])}.")
    print(f"So, the factorization of the sample space size is {' * '.join([f'{p}^{e * num_dice}' for p, e in side_factors.items()])}.")
    print("The maximum number of mutually independent events is the sum of the exponents:")
    print(f"{' + '.join(equation_parts)} = {max_m}")

solve_max_independent_events()
