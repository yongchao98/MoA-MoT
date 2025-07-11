import math

# The size of the semidihedral group is 512.
group_order = 512

# For a semidihedral group SD_{2^n}, the order is 2^n.
# So, 2^n = 512, which gives n=9.
n = int(math.log2(group_order))

# The group presentation is <r, s | r^(2^(n-1)) = s^2 = 1, ...>
# The order of the cyclic subgroup <r> is 2^(n-1).
r_order = 2**(n - 1)

# Case 1: k is an odd integer.
# For any odd k, the power subgroup G^k is the group G itself.
# This contributes 1 unique subgroup to our total count.
num_from_odd_k = 1

# Case 2: k is an even integer.
# For any even k, the power subgroup G^k is a cyclic subgroup <r^d>,
# where d = gcd(k, r_order). Since k is even, d must also be an
# even divisor of r_order. We count the number of such divisors.
num_from_even_k = 0
for i in range(1, r_order + 1):
  # Check if i is a divisor of r_order
  if r_order % i == 0:
    # Check if the divisor is even
    if i % 2 == 0:
      num_from_even_k += 1

# The total number of power subgroups is the sum from the two cases.
total_subgroups = num_from_odd_k + num_from_even_k

print("The semidihedral group of size 512 is SD_512.")
print("The analysis of its power subgroups leads to two cases:")
print(f"1. For odd k: This yields {num_from_odd_k} unique subgroup (the group itself).")
print(f"2. For even k: This yields {num_from_even_k} unique subgroups (the even divisors of {r_order}).")
print("\nThe final equation for the total number of power subgroups is:")
print(f"{num_from_odd_k} + {num_from_even_k} = {total_subgroups}")