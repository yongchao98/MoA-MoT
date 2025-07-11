import math

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors and their exponents for a given number n.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def main():
    """
    Calculates the largest possible number of mutually independent events
    for rolling 100 six-sided dice.
    """
    num_dice = 100
    num_sides = 6

    # The size of the sample space is num_sides ^ num_dice.
    # We need the prime factorization of this size.
    # Size = 6^100 = (2 * 3)^100 = 2^100 * 3^100.
    
    # We can get the exponents directly from the prime factors of the number of sides.
    side_factors = get_prime_factorization(num_sides)
    
    total_exponents = {}
    for prime, exponent in side_factors.items():
        total_exponents[prime] = exponent * num_dice

    # The maximum number of mutually independent events (m) is the sum of these exponents.
    m = sum(total_exponents.values())
    
    print("The experiment is rolling 100 dice, each with 6 sides.")
    print(f"The number of sides on a die is {num_sides}.")
    print(f"The prime factorization of {num_sides} is: {side_factors}")
    print("The size of the sample space is 6^100. Its prime factorization is 2^100 * 3^100.")
    print("The largest possible value of m is the sum of the exponents of the prime factors of the sample space size.")
    
    exponents = list(total_exponents.values())
    equation_parts = [str(e) for e in exponents]
    
    print("The exponents are:", ", ".join(equation_parts))
    print(f"The final equation is: m = {' + '.join(equation_parts)}")
    
    # We output each number in the final equation.
    # Since the exponents might change with input, this loop is more general.
    final_equation_str = "m = "
    for i, exp in enumerate(exponents):
        final_equation_str += str(exp)
        if i < len(exponents) - 1:
            final_equation_str += " + "
    final_equation_str += f" = {m}"
    print(final_equation_str)

main()
