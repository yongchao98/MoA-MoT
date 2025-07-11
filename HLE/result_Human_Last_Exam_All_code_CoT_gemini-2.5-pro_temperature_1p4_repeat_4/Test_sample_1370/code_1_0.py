import collections

def get_prime_factorization(n):
    """
    Returns a dictionary of the prime factorization of n.
    For example, get_prime_factorization(12) returns {2: 2, 3: 1}
    """
    factors = collections.defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def solve():
    """
    Calculates the largest possible number of mutually independent events
    for rolling a number of dice with a given number of sides.
    """
    num_dice = 100
    num_sides = 6

    # According to a theorem in combinatorics, the maximum number of mutually
    # independent events in an equiprobable sample space of size n is the sum
    # of the exponents in the prime factorization of n.
    # Here, n = num_sides ^ num_dice.

    # Step 1: Find prime factorization of the number of sides.
    side_factors = get_prime_factorization(num_sides)

    # Step 2: The exponents in the factorization of n = (s)^d are the exponents
    # of the factors of s, each multiplied by d. Sum these up.
    exponents = []
    for p, a in side_factors.items():
        exponent = a * num_dice
        exponents.append(exponent)
    
    m = sum(exponents)

    # Step 3: Print the result and the final equation.
    print(f"The number of sides on each die is {num_sides}.")
    print(f"The number of dice is {num_dice}.")
    print(f"The size of the sample space is n = {num_sides}^{num_dice}.")
    
    side_factor_str = " * ".join([f"{p}^{a}" for p, a in side_factors.items()])
    print(f"\nThe prime factorization of the number of sides ({num_sides}) is: {side_factor_str}.")
    
    n_factor_str = " * ".join([f"{p}^{a*num_dice}" for p, a in side_factors.items()])
    print(f"Therefore, the prime factorization of n is: {n_factor_str}.")
    
    print("\nThe largest possible number of mutually independent events (m) is the sum of these exponents.")
    
    equation_parts = [str(e) for e in exponents]
    equation_str = " + ".join(equation_parts)
    
    print(f"m = {equation_str} = {m}")

solve()