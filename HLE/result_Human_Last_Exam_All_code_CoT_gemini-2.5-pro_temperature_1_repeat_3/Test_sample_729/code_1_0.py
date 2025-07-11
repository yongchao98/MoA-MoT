import math

# The problem asks for the number of power subgroups in the
# generalized quaternion group of size 128.

# Let the group be G. The size of G is |G| = 128.
# The generalized quaternion group of order 2^n is denoted Q_{2^n}.
# To find n, we solve the equation 2^n = 128.
group_order = 128
n = int(math.log2(group_order))

# For a generalized quaternion group Q_{2^n} (with n >= 3),
# the number of distinct power subgroups is known to be exactly n.
# A power subgroup is a subgroup of the form G^k = {g^k for g in G}.
#
# The n distinct power subgroups for Q_{2^n} are:
# G^1, G^2, G^4, ..., G^(2^(n-1))
# where G^(2^(n-1)) is the trivial subgroup {1}.

# For G = Q_128, we have calculated n = 7.
number_of_power_subgroups = n

print(f"The generalized quaternion group of size {group_order} is of the form Q_(2^n).")
print(f"Solving 2^n = {group_order}, we find n = {n}.")
print(f"The number of power subgroups in Q_(2^n) is equal to n.")
print(f"Therefore, the number of power subgroups for the group of size {group_order} is {number_of_power_subgroups}.")
