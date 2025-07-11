import math

def get_prime_factorization_exponents(n):
    """
    Calculates the exponents of the prime factors of a number n.
    Returns a dictionary of {prime: exponent}.
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
    Calculates the maximum number of mutually independent events
    for rolling a set of dice.
    """
    num_dice = 100
    num_sides = 6

    # The theorem states that the max number of independent events is the
    # sum of the exponents in the prime factorization of the sample space size N.
    # N = num_sides ^ num_dice
    # Prime factorization of N can be found from the prime factorization of num_sides.
    
    # 1. Get prime factors of the number of sides.
    side_factors = get_prime_factorization_exponents(num_sides)
    
    # 2. The exponents of the prime factors of N = num_sides^num_dice are
    # the exponents of num_sides' factors multiplied by num_dice.
    total_exponents = []
    for prime in side_factors:
        exponent_in_N = side_factors[prime] * num_dice
        total_exponents.append(exponent_in_N)
        
    # 3. The maximum number of events 'm' is the sum of these exponents.
    max_m = sum(total_exponents)

    # 4. Print the final equation as requested.
    exponents_str = " + ".join(map(str, total_exponents))
    print(f"The size of the sample space is {num_sides}^{num_dice}.")
    print(f"The prime factorization of {num_sides} is {' * '.join([f'{p}^{e}' for p, e in side_factors.items()])}.")
    print(f"The prime factorization of the sample space size is {' * '.join([f'{p}^{e}' for p, e in zip(side_factors.keys(), total_exponents)])}.")
    print("\nThe maximum number of mutually independent events (m) is the sum of the exponents of this prime factorization.")
    print("\nThe final equation is:")
    print(f"m = {exponents_str} = {max_m}")

solve_max_independent_events()

<<<200>>>