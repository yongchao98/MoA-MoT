import math

def get_sum_of_prime_exponents(n):
    """
    Calculates the sum of the exponents in the prime factorization of a number n.
    For example, for n = 12 = 2^2 * 3^1, it returns 2 + 1 = 3.
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

def solve_dice_problem():
    """
    Solves the problem of finding the maximum number of mutually independent events
    for rolling a number of n-sided dice.
    """
    num_dice = 100
    num_sides = 6

    print("Step 1: The problem is to find the maximum number of mutually independent events.")
    print("This is equivalent to the sum of the exponents in the prime factorization of the sample space size, N.")
    
    print(f"\nStep 2: Calculate the size of the sample space, N.")
    print(f"N = (number of sides) ^ (number of dice) = {num_sides}^{num_dice}.")

    print(f"\nStep 3: Find the prime factorization of the number of sides, {num_sides}.")
    prime_factors_of_sides = get_sum_of_prime_exponents(num_sides)
    
    factor_str_parts = []
    exponents = []
    for base, exponent in prime_factors_of_sides.items():
        factor_str_parts.append(f"{base}^{exponent}")
        exponents.append(exponent)
    factor_str = " * ".join(factor_str_parts)
    print(f"The prime factorization of {num_sides} is {factor_str}.")

    sum_of_exponents_one_die = sum(exponents)
    print(f"The sum of the exponents for a single die is: {' + '.join(map(str, exponents))} = {sum_of_exponents_one_die}")
    
    print(f"\nStep 4: The prime factorization of N = {num_sides}^{num_dice} is ({factor_str})^{num_dice}.")
    print("The exponents in the prime factorization of N are each original exponent multiplied by the number of dice.")
    
    max_m = num_dice * sum_of_exponents_one_die
    
    print("\nStep 5: The largest possible value of m is the sum of these new exponents.")
    print("This is calculated by multiplying the number of dice by the sum of exponents for a single die.")
    
    # Printing the final equation with all numbers
    # Final equation format: num_dice * (sum of exponents for one die) = result
    print("\nFinal Calculation:")
    final_equation_str = f"{num_dice} * ({' + '.join(map(str, exponents))}) = {max_m}"
    print(final_equation_str)

solve_dice_problem()
<<<200>>>