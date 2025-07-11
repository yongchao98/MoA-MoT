import math

# Problem parameters
# D = (C_2)^5, so D is a vector space over F_2 of dimension 5.
dim_D = 5
p = 2
# Order of the defect group D
order_D = p**dim_D
# Order of the inertial quotient E
order_E = 5

# --- Step 1: Calculate l(B) ---
# l(B) is the number of irreducible Brauer characters.
# For a block with an abelian defect group, l(B) is the number of
# conjugacy classes of the inertial quotient E.
# Since E has prime order 5, it is a cyclic group C_5.
# The number of conjugacy classes in a cyclic group is its order.
l_B = order_E

# --- Step 2: Calculate k(B) ---
# k(B) is the number of irreducible ordinary characters.
# For a block with an abelian defect group, k(B) is the number of orbits
# of E acting on the character group of D, denoted D-hat.
# D-hat is isomorphic to D, so |D-hat| = 2^5 = 32.
size_D_hat = order_D

# We use Burnside's Lemma to count the orbits:
# k(B) = (1/|E|) * sum_{g in E} |Fix(g)|
# where Fix(g) is the set of elements in D-hat fixed by g.

# The identity element of E fixes all elements of D-hat.
num_fixed_by_identity = size_D_hat

# For any non-identity element g in E, we determine the number of fixed points.
# The action of E on D-hat is a 5-dimensional representation of C_5 over F_2.
# This representation decomposes into a 1-dim trivial and a 4-dim irreducible part.
# The fixed points for a non-identity element correspond to the trivial part.
# The dimension of the fixed-point subspace is 1.
dim_fixed_space = 1
num_fixed_by_non_identity = p**dim_fixed_space

# There are |E| - 1 non-identity elements in E.
num_non_identity_elements = order_E - 1

# Apply Burnside's Lemma
sum_of_fixed_points = num_fixed_by_identity + num_non_identity_elements * num_fixed_by_non_identity
k_B = sum_of_fixed_points / order_E
# We expect an integer result
k_B_int = int(k_B)

# --- Step 3: Compute and print the final result ---
result = k_B_int - l_B

print(f"Based on the theory of blocks with abelian defect groups:")
print(f"1. Calculation of l(B):")
print(f"l(B) = number of conjugacy classes of E. Since |E| = {order_E}, E is cyclic.")
print(f"So, l(B) = {l_B}.")
print("")
print(f"2. Calculation of k(B):")
print(f"k(B) = number of E-orbits on D-hat.")
print(f"Using Burnside's Lemma: k(B) = (1/{order_E}) * (|Fix(id)| + {num_non_identity_elements} * |Fix(g)|).")
print(f"   |Fix(id)| = |D-hat| = 2^{dim_D} = {size_D_hat}.")
print(f"   |Fix(g)| for g != id is 2^{dim_fixed_space} = {num_fixed_by_non_identity}.")
print(f"k(B) = ({num_fixed_by_identity} + {num_non_identity_elements} * {num_fixed_by_non_identity}) / {order_E} = {sum_of_fixed_points} / {order_E} = {k_B_int}")
print("")
print(f"3. Final result:")
print(f"k(B) - l(B) = {k_B_int} - {l_B} = {result}")
