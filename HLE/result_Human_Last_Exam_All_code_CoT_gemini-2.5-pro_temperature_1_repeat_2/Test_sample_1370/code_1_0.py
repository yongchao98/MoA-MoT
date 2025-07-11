import collections

def get_prime_factorization(n):
    """
    Returns a dictionary mapping each prime factor to its exponent.
    For example, get_prime_factorization(12) returns {2: 2, 3: 1}.
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

def solve_dice_independence():
    """
    Calculates the largest possible number of mutually independent events
    for rolling 100 6-sided dice.
    """
    num_dice = 100
    num_sides = 6

    # The maximum number of mutually independent events on an equiprobable sample space of size N
    # is the sum of the exponents in the prime factorization of N.
    # Here, N = num_sides ^ num_dice = 6^100.

    # 1. Find the prime factorization of the number of sides.
    base_factors = get_prime_factorization(num_sides)

    # 2. Calculate the exponents for the prime factorization of the total sample space size (6^100).
    final_exponents = []
    # Sorting keys to ensure a consistent order, e.g., 2^100 * 3^100
    for prime in sorted(base_factors.keys()):
        exponent = base_factors[prime]
        final_exponents.append(exponent * num_dice)

    # 3. The maximum number of events is the sum of these final exponents.
    max_m = sum(final_exponents)

    # 4. Print the final calculation as an equation.
    equation_str = " + ".join(map(str, final_exponents))
    print("The size of the sample space is 6^100.")
    print("The prime factorization of 6^100 is 2^100 * 3^100.")
    print("The largest possible value of m is the sum of these exponents.")
    print(f"m = {equation_str} = {max_m}")

solve_dice_independence()