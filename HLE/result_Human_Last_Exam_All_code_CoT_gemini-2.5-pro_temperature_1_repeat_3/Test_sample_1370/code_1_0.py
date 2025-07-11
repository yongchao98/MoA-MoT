import math

def get_prime_factorization(n):
    """
    Returns a dictionary representing the prime factorization of n.
    e.g., get_prime_factorization(12) returns {2: 2, 3: 1} for 12 = 2^2 * 3^1
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

def solve_dice_problem():
    """
    Calculates the largest possible number of mutually independent events
    for rolling 100 six-sided dice.
    """
    num_dice = 100
    num_sides = 6

    print(f"The experiment involves rolling {num_dice} dice, each with {num_sides} sides.")
    print(f"The size of the sample space is {num_sides}^{num_dice}.")
    print("-" * 20)

    # Step 1: Find the prime factorization of the number of sides.
    side_factors = get_prime_factorization(num_sides)
    
    print(f"The prime factorization of the number of sides ({num_sides}) is:")
    factor_str = " * ".join([f"{p}^{e}" for p, e in side_factors.items()])
    print(f"{num_sides} = {factor_str}")
    print("-" * 20)

    # Step 2: Determine the exponents for the prime factorization of the total sample space size.
    # The exponents of the total space size are the exponents of the side factors multiplied by the number of dice.
    total_space_exponents = {p: e * num_dice for p, e in side_factors.items()}
    
    print(f"The prime factorization of the total sample space size ({num_sides}^{num_dice}) is:")
    total_factor_str = " * ".join([f"{p}^{e}" for p, e in total_space_exponents.items()])
    print(f"{num_sides}^{num_dice} = {total_factor_str}")
    print("-" * 20)

    # Step 3: Sum the exponents to find the maximum number of independent events.
    # This is based on the theorem for finite equiprobable spaces.
    max_m = sum(total_space_exponents.values())
    
    exponents = list(total_space_exponents.values())
    
    print("According to the theorem, the maximum number of mutually independent events")
    print("is the sum of the exponents in the prime factorization of the sample space size.")
    print("\nThe final calculation is:")
    
    # Building the final equation string
    equation_parts = [str(e) for e in exponents]
    equation_str = " + ".join(equation_parts)
    
    print(f"{equation_str} = {max_m}")
    
    print("\nTherefore, the largest possible value of m is:")
    print(max_m)

solve_dice_problem()