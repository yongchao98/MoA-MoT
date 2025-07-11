def solve_dice_problem():
    """
    Calculates the largest possible number of mutually independent events
    when rolling 100 6-sided dice.
    """
    num_sides = 6
    num_dice = 100

    def get_prime_factorization(n):
        """Returns the prime factorization of n as a dictionary {prime: exponent}."""
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

    # The size of the sample space is |Ω| = num_sides ^ num_dice.
    # A theorem in probability states that the maximum number of mutually independent
    # events is the sum of the exponents in the prime factorization of |Ω|.
    #
    # First, we find the prime factorization of the base (num_sides).
    # For num_sides = 6, this is {2: 1, 3: 1}.
    base_factors = get_prime_factorization(num_sides)

    # For |Ω| = 6^100, the prime factorization is 2^100 * 3^100.
    # The exponents are derived by multiplying the base exponents by num_dice.
    final_exponents = []
    for prime in sorted(base_factors.keys()):
        base_exponent = base_factors[prime]
        final_exponent = base_exponent * num_dice
        final_exponents.append(final_exponent)

    # The maximum number of events 'm' is the sum of these final exponents.
    m = sum(final_exponents)

    # Print the final calculation as requested.
    equation_parts = [str(e) for e in final_exponents]
    equation_str = " + ".join(equation_parts)
    print(f"The maximum number of events is the sum of the exponents in the prime factorization of the sample space size (6^100 = 2^100 * 3^100).")
    print(f"The final equation is: {equation_str} = {m}")

solve_dice_problem()