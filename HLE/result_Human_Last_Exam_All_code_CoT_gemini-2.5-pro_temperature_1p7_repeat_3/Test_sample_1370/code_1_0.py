import collections

def prime_factorize(n):
    """
    Returns a dictionary of prime factors of n and their exponents.
    For example, prime_factorize(12) returns {2: 2, 3: 1}.
    """
    factors = collections.defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def solve_dice_problem():
    """
    Solves the problem of finding the largest possible number of mutually
    independent events for rolling 100 six-sided dice.
    """
    num_dice = 100
    sides_per_die = 6

    # The size of the sample space is N = sides_per_die ** num_dice.
    # The maximum number of mutually independent events is the sum of the
    # exponents in the prime factorization of N.
    # This can be found by prime factorizing the base (sides_per_die) and
    # multiplying the exponents by the number of dice.

    base = sides_per_die
    exponent = num_dice

    # Get prime factorization of the base, e.g., 6 -> {2: 1, 3: 1}
    base_factors = prime_factorize(base)
    
    total_exponents = 0
    exponent_list = []
    
    # Calculate the final exponents for N = 6^100 and add them to a list.
    # For each prime p with exponent 'exp' in the base factorization,
    # its exponent in the final number is 'exp * exponent'.
    for p, exp in sorted(base_factors.items()):
        final_exponent = exp * exponent
        exponent_list.append(final_exponent)
        total_exponents += final_exponent
        
    # Format the final equation for printing.
    equation_str = " + ".join(map(str, exponent_list))
    
    print("The largest possible value of m is the sum of the exponents in the prime factorization of the sample space size (6^100).")
    print("The final equation is:")
    print(f"{equation_str} = {total_exponents}")

solve_dice_problem()