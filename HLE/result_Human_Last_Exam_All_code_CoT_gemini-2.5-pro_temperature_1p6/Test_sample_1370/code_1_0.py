def solve_dice_problem():
    """
    Calculates the largest possible number of mutually independent events
    when rolling 100 regular 6-sided dice.
    """

    # The number of dice rolled.
    num_dice = 100
    # The number of sides on each die.
    num_sides = 6

    # The size of the sample space is N = num_sides ^ num_dice.
    # We need to find the sum of the exponents in the prime factorization of N.
    # This is equivalent to num_dice * (sum of exponents in prime factorization of num_sides).

    # Prime factorization of the number of sides (6).
    # 6 = 2 * 3. The exponents are 1 for prime 2, and 1 for prime 3.
    exponent_of_2_in_base = 1
    exponent_of_3_in_base = 1

    # For N = 6^100, the exponents are multiplied by 100.
    # N = (2^1 * 3^1)^100 = 2^100 * 3^100.
    exponent_of_2_in_N = exponent_of_2_in_base * num_dice
    exponent_of_3_in_N = exponent_of_3_in_base * num_dice

    # The maximum number of mutually independent events is the sum of these exponents.
    max_m = exponent_of_2_in_N + exponent_of_3_in_N

    print("The largest possible value of m is the sum of the exponents in the prime factorization of the size of the sample space, N = 6^100.")
    print("The prime factorization of N is 2^100 * 3^100.")
    print("The final calculation is the sum of the exponents:")
    print(f"{exponent_of_2_in_N} + {exponent_of_3_in_N} = {max_m}")

solve_dice_problem()