import math

def get_prime_factorization_exponents(n):
    """
    Returns a dictionary of prime factors and their exponents.
    e.g., get_prime_factorization_exponents(12) returns {2: 2, 3: 1}
    """
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

# Problem parameters from the user's request
num_dice = 100
num_sides = 6

print(f"We are rolling {num_dice} dice, each with {num_sides} sides.")
print(f"The size of the sample space is {num_sides}^{num_dice}.")

# The maximum number of mutually independent events is the sum of the exponents
# in the prime factorization of the size of the sample space.
# N = (num_sides)^num_dice = (p1^a1 * p2^a2 * ...)^num_dice = p1^(a1*num_dice) * p2^(a2*num_dice) * ...
# m = (a1 * num_dice) + (a2 * num_dice) + ...

# 1. Find the prime factorization of the number of sides.
prime_exponents_of_sides = get_prime_factorization_exponents(num_sides)

print(f"The prime factorization of the number of sides ({num_sides}) gives us the base exponents:")
for p, exp in prime_exponents_of_sides.items():
    print(f"  - Prime: {p}, Exponent: {exp}")

# 2. Calculate the maximum number of independent events.
max_m = 0
calculation_parts = []
for prime, exponent in prime_exponents_of_sides.items():
    term = exponent * num_dice
    max_m += term
    # This part prepares the string for the final equation printout
    calculation_parts.append(f"{exponent} * {num_dice}")

equation_str = " + ".join(calculation_parts)

print("\nThe largest possible value of m is the sum of these exponents, each multiplied by the number of dice.")
print(f"The final equation is:")
print(f"m = {equation_str}")
print(f"m = {max_m}")