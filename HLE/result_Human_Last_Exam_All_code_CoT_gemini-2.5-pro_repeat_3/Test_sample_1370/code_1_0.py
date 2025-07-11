import collections

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors and their exponents for a given number n.
    """
    factors = collections.defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
       factors[temp] += 1
    return factors

def solve_max_independent_events():
    """
    Calculates the maximum number of mutually independent events for rolling
    a number of dice with a given number of sides.
    """
    num_dice = 100
    num_sides = 6

    print("This script calculates the largest possible number of mutually independent events (m).")
    print("---")

    # Step 1: Define the size of the sample space
    print("Step 1: Determine the size of the sample space (N).")
    print(f"The experiment involves rolling {num_dice} dice, each with {num_sides} sides.")
    print(f"Total outcomes N = {num_sides}^{num_dice}.")
    print("---")

    # Step 2: Find the prime factorization of N
    print("Step 2: Find the prime factorization of N.")
    print("According to a theorem in combinatorics, m is the sum of the exponents in the prime factorization of N.")
    
    base_factors = get_prime_factorization(num_sides)
    
    print(f"The prime factorization of the base ({num_sides}) is: {' * '.join([f'{p}^{e}' for p, e in base_factors.items()])}")
    
    exponents_in_N = {p: e * num_dice for p, e in base_factors.items()}
    factor_strings = [f"{p}^{e}" for p, e in exponents_in_N.items()]
    
    print(f"So, the prime factorization of N = {num_sides}^{num_dice} is: {' * '.join(factor_strings)}.")
    print("---")

    # Step 3: Sum the exponents
    print("Step 3: Sum the exponents to find m.")
    exponent_values = list(exponents_in_N.values())
    max_m = sum(exponent_values)
    
    equation_parts = [str(e) for e in exponent_values]
    print(f"The exponents from the factorization of N are: {', '.join(equation_parts)}.")
    print("The final equation for m is the sum of these exponents:")
    print(f"m = {' + '.join(equation_parts)} = {max_m}")
    print("---")
    
    print(f"The largest possible value of m is {max_m}.")

solve_max_independent_events()