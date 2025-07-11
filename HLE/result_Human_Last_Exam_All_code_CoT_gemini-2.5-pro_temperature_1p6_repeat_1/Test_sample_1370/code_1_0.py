import math

def get_prime_factorization(n):
    """
    Returns a dictionary where keys are prime factors and values are their exponents.
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

def solve():
    """
    Calculates the largest possible number of mutually independent events
    for rolling 100 6-sided dice.
    """
    num_dice = 100
    num_sides = 6

    # The size of the sample space is N = num_sides ^ num_dice, which is 6^100.
    # A theorem states that the max number of mutually independent events is the
    # sum of the exponents in the prime factorization of the sample space size.

    # 1. Get the prime factorization of the number of sides.
    side_factors = get_prime_factorization(num_sides)
    
    # 2. The exponents in the prime factorization of N = (num_sides)^num_dice are
    #    the exponents from side_factors multiplied by num_dice.
    total_exponents = []
    for prime, exponent in side_factors.items():
        total_exponents.append(exponent * num_dice)

    # 3. Sum these exponents to get the maximum number of events, m.
    m = sum(total_exponents)

    print("The problem is to find the maximum number of mutually independent events (m).")
    print("This number is the sum of the exponents of the prime factorization of the sample space size, |Ω|.")
    print(f"For {num_dice} dice with {num_sides} sides, the size is |Ω| = {num_sides}^{num_dice}.")
    print(f"The prime factorization of {num_sides} is {side_factors}.")
    print(f"So, the prime factorization of |Ω| is ({' * '.join(map(str, side_factors.keys()))})^{num_dice}, which gives exponents: {total_exponents}")
    
    # 4. Print the final calculation step-by-step.
    equation_parts = [str(exp) for exp in total_exponents]
    final_equation = f"m = {' + '.join(equation_parts)} = {m}"
    print("\nThe final calculation is:")
    print(final_equation)

solve()
<<<200>>>