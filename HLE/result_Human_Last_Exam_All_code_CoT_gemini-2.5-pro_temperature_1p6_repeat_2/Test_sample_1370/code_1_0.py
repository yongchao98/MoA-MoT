import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary of {prime: exponent}.
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
    Calculates the largest possible number of mutually independent events
    for rolling a set of dice and prints the step-by-step logic.
    """
    num_dice = 100
    num_sides = 6

    print(f"This problem involves {num_dice} dice, each with {num_sides} sides.")
    print("The total size of the sample space is N = (number of sides)^(number of dice).")
    print("A theorem states that the max number of mutually independent events is the sum of the exponents in the prime factorization of N.")
    print("-" * 40)

    # Step 1: Find prime factors of the number of sides.
    print("Step 1: Find the prime factorization of the number of sides on a single die.")
    side_factors = get_prime_factorization(num_sides)
    
    side_factor_parts = []
    for prime, exponent in side_factors.items():
        if exponent > 1:
            side_factor_parts.append(f"{prime}^{exponent}")
        else:
            side_factor_parts.append(str(prime))
            
    print(f"The number of sides is {num_sides}.")
    print(f"The prime factorization of {num_sides} is: {' * '.join(side_factor_parts)}")
    print("-" * 40)

    # Step 2: Determine exponents for the full sample space size N.
    print("Step 2: Determine the exponents in the prime factorization of the total sample space size, N.")
    print(f"N = {num_sides}^{num_dice} = ({' * '.join(side_factor_parts)})^{num_dice}")
    
    final_exponents = []
    final_factor_parts = []
    for prime, exponent in side_factors.items():
        final_exponent = exponent * num_dice
        final_exponents.append(final_exponent)
        final_factor_parts.append(f"{prime}^{final_exponent}")
    
    print(f"The prime factorization of N is: {' * '.join(final_factor_parts)}")
    print(f"The exponents of the prime factors are: {', '.join(map(str, final_exponents))}")
    print("-" * 40)

    # Step 3: Sum the exponents to get the final answer.
    print("Step 3: The largest possible value of m is the sum of these exponents.")
    m = sum(final_exponents)
    
    equation_parts = [str(exp) for exp in final_exponents]
    print(f"The final equation is: m = {' + '.join(equation_parts)}")
    print(f"The largest possible value of m is: {m}")

solve_max_independent_events()