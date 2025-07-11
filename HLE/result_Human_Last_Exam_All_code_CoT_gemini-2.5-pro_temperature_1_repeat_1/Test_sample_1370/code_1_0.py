def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary of {prime: exponent}.
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

def solve_dice_independence():
    """
    Calculates the largest possible number of mutually independent events
    for rolling a given number of dice with a given number of sides.
    """
    num_dice = 100
    num_sides = 6

    # The size of the sample space is |Ω| = num_sides ^ num_dice.
    # The maximum number of mutually independent events is the sum of the exponents
    # in the prime factorization of |Ω|. This can be found by finding the prime
    # factorization of num_sides and multiplying each exponent by num_dice.

    print(f"For {num_dice} dice with {num_sides} sides each:")
    
    # 1. Find prime factorization of the number of sides.
    side_factors = get_prime_factorization(num_sides)
    
    # 2. Calculate the exponents for the total sample space |Ω| and sum them up.
    m_terms = []
    max_m = 0
    
    # Sort the primes for a consistent output order.
    sorted_primes = sorted(side_factors.keys())

    for prime in sorted_primes:
        exponent = side_factors[prime]
        term = exponent * num_dice
        m_terms.append(term)
        max_m += term
    
    # 3. Print the step-by-step reasoning and the final result.
    print(f"The size of the sample space is |Ω| = {num_sides}^{num_dice}.")
    prime_factor_str = " * ".join([f"{p}^{e}" for p, e in side_factors.items()])
    print(f"The prime factorization of {num_sides} is {prime_factor_str}.")
    exponent_str = " * ".join([f"{p}^{e * num_dice}" for p, e in side_factors.items()])
    print(f"Thus, the prime factorization of |Ω| is {exponent_str}.")
    print("The maximum number of mutually independent events is the sum of these exponents.")
    
    equation_str = " + ".join(map(str, m_terms))
    print("\nThe largest possible value of m is:")
    print(f"{equation_str} = {max_m}")

solve_dice_independence()