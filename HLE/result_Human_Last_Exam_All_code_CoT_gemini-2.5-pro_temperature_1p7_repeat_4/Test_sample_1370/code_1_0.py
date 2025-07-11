import math

def get_sum_of_prime_exponents(n):
    """
    Calculates the sum of the exponents in the prime factorization of n.
    For n = p1^a1 * p2^a2 * ..., returns the sum a1 + a2 + ...
    Also returns string representations for the factorization and the sum of exponents.
    """
    original_n = n
    sum_exponents = 0
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            sum_exponents += 1
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        sum_exponents += 1
        factors[temp_n] = factors.get(temp_n, 0) + 1
    
    exponent_sum_str = " + ".join([str(a) for p, a in sorted(factors.items())])
    
    return sum_exponents, exponent_sum_str

def solve_dice_problem():
    """
    Calculates the largest possible number of mutually independent events
    for rolling k n-sided dice.
    """
    # Number of dice
    num_dice = 100

    # Number of sides on each die
    num_sides = 6

    # Calculate the sum of exponents for the number of sides
    sum_exponents, exponent_sum_str = get_sum_of_prime_exponents(num_sides)

    # The maximum number of mutually independent events m is num_dice * sum_exponents
    m_max = num_dice * sum_exponents

    # Print the explanation and the result as a full equation
    print("The largest possible value of m is found by the formula:")
    print("m = (Number of dice) * (Sum of exponents in the prime factorization of the number of sides)")
    print("\nCalculating the components:")
    print(f"Number of dice = {num_dice}")
    print(f"Number of sides = {num_sides}")
    print(f"The sum of exponents for the prime factors of {num_sides} (which is 2^1 * 3^1) is {exponent_sum_str} = {sum_exponents}.")
    print("\nThe final equation with all numbers is:")
    print(f"{m_max} = {num_dice} * {sum_exponents}")

solve_dice_problem()