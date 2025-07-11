import math

def get_prime_factorization_with_exponents(n):
    """
    Calculates the prime factorization of a number n.
    Returns a dictionary where keys are prime factors and values are their exponents.
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
    Calculates and prints the largest possible number of mutually independent events
    for rolling a number of dice.
    """
    num_dice = 100
    sides_per_die = 6

    print("This problem can be solved by finding the total number of prime factors (with multiplicity) in the size of the sample space.")
    print(f"The sample space consists of all outcomes from rolling {num_dice} {sides_per_die}-sided dice.")
    print(f"The total number of outcomes is |Ω| = {sides_per_die}^{num_dice}.")
    print("\nStep 1: Find the prime factorization of the number of sides of a single die.")
    
    side_factors = get_prime_factorization_with_exponents(sides_per_die)
    
    factor_str_parts = []
    for prime, exponent in side_factors.items():
        factor_str_parts.append(f"{prime}^{exponent}")
    
    print(f"The prime factorization of {sides_per_die} is: {' * '.join(factor_str_parts)}")

    print("\nStep 2: The number of prime factors of |Ω| is the number of dice multiplied by the sum of the exponents of the prime factors of a single die's sides.")

    sum_of_exponents = sum(side_factors.values())
    exponent_sum_str = " + ".join(map(str, side_factors.values()))
    
    max_m = num_dice * sum_of_exponents
    
    print("\nStep 3: Calculate the final value 'm'.")
    print(f"m = (number of dice) * (sum of exponents of prime factors of {sides_per_die})")
    print(f"m = {num_dice} * ({exponent_sum_str})")
    print(f"m = {num_dice} * {sum_of_exponents}")
    print(f"m = {max_m}")
    
    print("\nThus, the largest possible value of m is:")
    print(max_m)

solve_max_independent_events()
