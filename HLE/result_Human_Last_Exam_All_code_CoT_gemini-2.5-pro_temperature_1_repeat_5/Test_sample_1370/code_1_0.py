def get_prime_factorization_exponents(n):
    """
    Finds the exponents of the prime factors of a number n.
    For example, for n = 12 = 2^2 * 3^1, it returns [2, 1].
    """
    factors = []
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            count = 0
            while temp_n % d == 0:
                count += 1
                temp_n //= d
            factors.append(count)
        d += 1
    if temp_n > 1:
        factors.append(1)
    return factors

def solve_dice_problem():
    """
    Calculates the largest possible number of mutually independent events
    for rolling a number of dice with a given number of sides.
    """
    num_dice = 100
    num_sides = 6

    # The total number of outcomes is N = num_sides ^ num_dice.
    # The max number of independent events is the sum of exponents in the prime factorization of N.
    # This is equivalent to finding the exponents for num_sides and multiplying by num_dice.

    # 1. Find the exponents of the prime factors of the base (num_sides).
    # For 6 = 2^1 * 3^1, the exponents are [1, 1].
    base_exponents = get_prime_factorization_exponents(num_sides)

    # 2. The exponents for the total number of outcomes (num_sides^num_dice)
    # are the base exponents each multiplied by num_dice.
    final_exponents = [e * num_dice for e in base_exponents]

    # 3. The maximum number of events is the sum of these final exponents.
    max_m = sum(final_exponents)

    # 4. Print the final equation and the result.
    equation_parts = [str(e) for e in final_exponents]
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {max_m}")

solve_dice_problem()