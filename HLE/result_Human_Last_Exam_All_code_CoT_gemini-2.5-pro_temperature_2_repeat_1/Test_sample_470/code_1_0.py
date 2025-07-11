import math

# Step 1: Define the parameters from the problem statement.
# The defect group D is (C_2)^5.
dim_V = 5
order_of_D = 2**dim_V

# The inertial quotient E has order 5. E is isomorphic to the cyclic group C_5.
order_of_E = 5

print("Problem Parameters:")
print(f"Order of the defect group D, |D| = 2^{dim_V} = {order_of_D}")
print(f"Order of the inertial quotient E, |E| = {order_of_E}\n")

# Step 2: Calculate l(B), the number of irreducible Brauer characters.
# For a block with an abelian defect group, l(B) = |E|.
l_B = order_of_E
print("Calculation of l(B):")
print(f"The number of irreducible Brauer characters, l(B), is the order of the inertial quotient.")
print(f"l(B) = |E| = {l_B}\n")

# Step 3: Calculate k(B), the number of irreducible ordinary characters.
# k(B) is the number of orbits of E on Irr(D). We use Burnside's Lemma.
# |Irr(D)| is isomorphic to D, so it has |D| = 32 elements.
num_elements = order_of_D

# Number of elements fixed by the identity element of E.
# The identity fixes all elements of Irr(D).
fixed_points_identity = num_elements

# Number of elements fixed by any non-identity element of E.
# As explained in the plan, for a non-identity element g in E (order 5), its action
# on V = (F_2)^5 has a 1-dimensional fixed-point subspace.
# This subspace contains 2^1 = 2 elements.
dim_fixed_space_non_identity = 1
fixed_points_non_identity = 2**dim_fixed_space_non_identity

# There is 1 identity element and (order_of_E - 1) non-identity elements in E.
num_non_identity_elements = order_of_E - 1

# Apply Burnside's Lemma: k(B) = (1/|E|) * sum over g in E of |Fix(g)|
sum_of_fixed_points = fixed_points_identity * 1 + fixed_points_non_identity * num_non_identity_elements
k_B = int(sum_of_fixed_points / order_of_E)

print("Calculation of k(B):")
print(f"The number of ordinary characters, k(B), is the number of orbits of E on Irr(D).")
print(f"Using Burnside's Lemma: k(B) = (1/|E|) * ( |Fix(1)| + sum_{{g!=1}} |Fix(g)| )")
print(f"|Fix(1)| (elements fixed by identity) = {fixed_points_identity}")
print(f"|Fix(g)| (elements fixed by a non-identity element) = {fixed_points_non_identity}")
print(f"Number of non-identity elements in E = {num_non_identity_elements}")
print(f"k(B) = (1/{order_of_E}) * ({fixed_points_identity} + {num_non_identity_elements} * {fixed_points_non_identity}) = (1/{order_of_E}) * ({sum_of_fixed_points}) = {k_B}\n")


# Step 4: Compute the final result k(B) - l(B).
result = k_B - l_B

print("Final Calculation:")
print(f"k(B) - l(B) = {k_B} - {l_B} = {result}")
