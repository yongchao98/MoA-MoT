def solve_max_independent_events():
    """
    Calculates the largest possible number of mutually independent events
    for rolling 100 6-sided dice.
    """
    num_dice = 100
    sides_per_die = 6

    # Prime factorization of the number of sides (6 = 2^1 * 3^1)
    base_exponents = {
        '2': 1,
        '3': 1
    }

    # The size of the sample space is 6^100.
    # The prime factorization of the sample space size is (2*3)^100 = 2^100 * 3^100.
    # The exponents in the prime factorization of the sample space size are:
    sample_space_exponents = {p: e * num_dice for p, e in base_exponents.items()}

    # The maximum number of mutually independent events is the sum of these exponents.
    exponent1 = sample_space_exponents['2']
    exponent2 = sample_space_exponents['3']
    max_m = exponent1 + exponent2
    
    # Example construction:
    # For each of the 100 dice, we can define two independent events:
    # 1. The outcome is even (probability 3/6 = 1/2).
    # 2. The outcome is a multiple of 3 (probability 2/6 = 1/3).
    # This gives 100 * 2 = 200 mutually independent events in total.
    
    print("The size of the sample space is 6^100, which has the prime factorization 2^100 * 3^100.")
    print("According to the theorem, the maximum number of mutually independent events is the sum of the exponents in the prime factorization of the sample space size.")
    print(f"The equation is: m = {exponent1} + {exponent2}")
    print(f"The largest possible value of m is: {max_m}")

solve_max_independent_events()
