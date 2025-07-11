import math

def get_prime_factorization_exponents(n):
    """
    Computes the exponents of the prime factors of a given number n.
    Returns a dictionary where keys are prime factors and values are their exponents.
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

# 1. Define the parameters of the problem
num_dice = 100
num_sides = 6

# 2. State the size of the sample space, N
# N is num_sides raised to the power of num_dice, i.e., 6^100.
# We don't need to compute this large number directly.
# We need its prime factorization.

print(f"The sample space consists of {num_sides}^{num_dice} equiprobable outcomes.")
print("We need to find the prime factorization of this number.")

# 3. Find the prime factorization of the base (number of sides)
base_factors = get_prime_factorization_exponents(num_sides)
print(f"The prime factorization of the number of sides ({num_sides}) is: {base_factors}")

# 4. Determine the exponents in the prime factorization of N = (base)^exponent
final_exponents = {}
for prime, exponent in base_factors.items():
    final_exponents[prime] = exponent * num_dice

print(f"The exponents in the prime factorization of N = {num_sides}^{num_dice} are: {final_exponents}")

# 5. Apply the theorem
# The maximum number of mutually independent events is the sum of these exponents.
max_m = sum(final_exponents.values())

print("\nA known theorem states that the maximum number of mutually independent events in an")
print("equiprobable sample space is the sum of the exponents of the prime factors of the")
print("sample space size.")

# 6. Display the final calculation clearly
equation_parts = [str(exp) for exp in final_exponents.values()]
equation_str = " + ".join(equation_parts)
print(f"\nThe largest possible value of m is the sum of these exponents:")
print(f"m = {equation_str}")
print(f"m = {max_m}")
