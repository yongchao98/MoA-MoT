import math

def euler_totient(n):
    """
    Calculates Euler's totient function, phi(n), which counts the positive
    integers up to a given integer n that are relatively prime to n.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    
    result = n
    p = 2
    temp_n = n
    while p * p <= temp_n:
        if temp_n % p == 0:
            while temp_n % p == 0:
                temp_n //= p
            result -= result // p
        p += 1
    if temp_n > 1:
        result -= result // temp_n
    return result

# The degree of the polynomial f(x) determines the order of its cyclic Galois group.
degree = 7
group_order = degree

# The number of generators of the cyclic group C_n is phi(n).
num_generators = euler_totient(group_order)

# The density of primes p for which f(x) remains irreducible mod p
# is the ratio of the number of generators to the group order.
density_numerator = num_generators
density_denominator = group_order

print(f"The polynomial is f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22.")
print(f"The Galois group of this polynomial is the cyclic group C_{degree}.")
print(f"The order of the Galois group is {group_order}.")
print(f"For f(x) to be irreducible modulo a prime p, the Frobenius element Frob_p must be a generator of the Galois group.")
print(f"The number of generators of C_{group_order} is given by Euler's totient function, phi({group_order}).")
print(f"phi({group_order}) = {num_generators}.")
print(f"The density is the ratio of the number of generators to the order of the group.")
print(f"Density = {density_numerator} / {density_denominator}")
result_fraction = f"{density_numerator}/{density_denominator}"
print(f"The final answer is {result_fraction}.")
