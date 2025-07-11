import math

# Group properties
group_order = 128
# Q_128 is the generalized quaternion group Q_{4n} with 4n = 128, so n = 32.
n_param = 32
# The order of the element x in the standard presentation is 2n.
order_x = 2 * n_param

print(f"The generalized quaternion group of size {group_order}, denoted Q_{group_order}, has the presentation:")
print(f"G = <x, y | x^{order_x} = 1, x^{n_param} = y^2, yx = x^-1 * y>\n")
print("A power subgroup H is a subgroup that can be written as H = G^k = {g^k | g in G} for some integer k.\n")

print("The distinct power subgroups are found by analyzing the set G^k for different types of k:")

# Case 1: k is odd
count_odd_k = 1
print(f"1. For any odd integer k, the set G^k is the group G itself. This gives {count_odd_k} unique subgroup: Q_128.")

# Case 2: k is congruent to 2 modulo 4
count_k_2_mod_4 = 1
print(f"2. For any integer k such that k = 2 (mod 4), G^k is the cyclic subgroup <x^2>. This gives {count_k_2_mod_4} unique subgroup: <x^2>.")

# Case 3: k is a multiple of 4
# For this case, the power subgroups are of the form <x^d>, where d must be a divisor of 64 and a multiple of 4.
# We will find these values of d.

# Find all divisors of order_x (64)
divisors_of_order_x = [i for i in range(1, order_x + 1) if order_x % i == 0]
# Filter for those that are multiples of 4
d_values = [d for d in divisors_of_order_x if d % 4 == 0]
count_k_4_mod_4 = len(d_values)

print(f"3. For any integer k that is a multiple of 4, G^k results in a cyclic subgroup of the form <x^d>, where d = gcd(k, {order_x}).")
print(f"   The possible values for d are the divisors of {order_x} that are also multiples of 4: {d_values}.")
print(f"   This gives {count_k_4_mod_4} distinct subgroups: <x^4>, <x^8>, <x^16>, <x^32>, and <x^64> (which is the trivial subgroup {{1}}).")

# Calculate the total count
total_count = count_odd_k + count_k_2_mod_4 + count_k_4_mod_4

print("\nTo find the total number of power subgroups, we sum the counts of unique subgroups from each case.")
print(f"The final calculation is: {count_odd_k} + {count_k_2_mod_4} + {count_k_4_mod_4} = {total_count}")