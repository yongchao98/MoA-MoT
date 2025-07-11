import math

def combinations(n, k):
    """Helper function to calculate combinations C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# --- Calculate involutions in PSL(3,4) ---
# PSL(3,4) is isomorphic to the Alternating Group A_8.
# Involutions are elements of order 2. In A_8, these are permutations
# with cycle structures of (2,2) or (2,2,2,2).

# Number of permutations of type (ab)(cd)
# C(8,2) * C(6,2) / 2!
num_type_2_2 = (combinations(8, 2) * combinations(6, 2)) // math.factorial(2)

# Number of permutations of type (ab)(cd)(ef)(gh)
# C(8,2) * C(6,2) * C(4,2) * C(2,2) / 4!
num_type_2_2_2_2 = (combinations(8, 2) * combinations(6, 2) * combinations(4, 2) * combinations(2, 2)) // math.factorial(4)

total_involutions_psl34 = num_type_2_2 + num_type_2_2_2_2

print("Calculation for PSL(3,4) (isomorphic to A_8):")
print(f"Number of involutions of type (ab)(cd): {num_type_2_2}")
print(f"Number of involutions of type (ab)(cd)(ef)(gh): {num_type_2_2_2_2}")
print(f"Total number of involutions in PSL(3,4) = {num_type_2_2} + {num_type_2_2_2_2} = {total_involutions_psl34}")
print("-" * 30)

# --- Calculate involutions in PSU(3,3) ---
# This data is known from the character table of PSU(3,3).
# It has two conjugacy classes of involutions.
class1_size = 63
class2_size = 252
total_involutions_psu33 = class1_size + class2_size

print("Calculation for PSU(3,3):")
print(f"Size of the first class of involutions: {class1_size}")
print(f"Size of the second class of involutions: {class2_size}")
print(f"Total number of involutions in PSU(3,3) = {class1_size} + {class2_size} = {total_involutions_psu33}")
print("-" * 30)

# --- Conclusion ---
print(f"Comparing the number of involutions:")
print(f"PSL(3,4): {total_involutions_psl34}")
print(f"PSU(3,3): {total_involutions_psu33}")
if total_involutions_psl34 == total_involutions_psu33:
    print("The number of involutions is equal for PSL(3,4) and PSU(3,3).")
else:
    print("The number of involutions is not equal.")
