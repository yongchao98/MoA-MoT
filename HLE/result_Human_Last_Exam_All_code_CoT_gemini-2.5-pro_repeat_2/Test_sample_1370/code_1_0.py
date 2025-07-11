def solve_dice_problem():
    """
    Calculates the largest possible number of mutually independent events
    for rolling 100 six-sided dice.
    """
    # Number of dice rolled
    num_dice = 100

    # The number of sides on each die
    num_sides = 6

    # The size of the sample space is num_sides ^ num_dice.
    # We need the prime factorization of this number.
    # Let's start with the prime factorization of the base (number of sides).
    # 6 = 2 * 3 = 2^1 * 3^1
    # The exponents in the prime factorization of the base are 1 and 1.
    base_exponents = [1, 1] 

    # For a sample space of size N = (p1^a1 * p2^a2 * ...)^k = p1^(k*a1) * p2^(k*a2) * ...,
    # the sum of exponents is k * (a1 + a2 + ...).
    # In our case, N = 6^100, so the exponents are 100 * 1 and 100 * 1.
    final_exponents = [exponent * num_dice for exponent in base_exponents]

    # The maximum number of mutually independent events is the sum of these final exponents.
    max_m = sum(final_exponents)
    
    # According to the problem instruction, we need to output each number in the final equation.
    print("The size of the sample space is 6^100.")
    print("The prime factorization of 6^100 is 2^100 * 3^100.")
    print("The exponents from the prime factorization are {} and {}.".format(final_exponents[0], final_exponents[1]))
    print("The largest possible value of m is the sum of these exponents:")
    print(f"{final_exponents[0]} + {final_exponents[1]} = {max_m}")

solve_dice_problem()