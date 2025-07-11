import math

# Step 1: Define the given parameters.
# D = (C_2)^5, so D is an elementary abelian group of order 2^5.
order_D = 2**5

# The inertial quotient E has order 5.
order_E = 5

print(f"The defect group D has order |D| = {order_D}.")
print(f"The inertial quotient E has order |E| = {order_E}.")
print("-" * 20)

# Step 2: Calculate l(B), the number of irreducible Brauer characters.
# For a block with an abelian defect group D and inertial quotient E,
# l(B) is the number of simple modules for the group algebra FE.
# Since char(F) = 2 does not divide |E| = 5, and F is large enough,
# this number is equal to the order of E.
l_B = order_E
print(f"The number of irreducible Brauer characters is l(B) = |E| = {l_B}.")
print("-" * 20)

# Step 3: Calculate k(B), the number of irreducible characters.
# For a block with an abelian defect group D, k(B) is the number of
# orbits of E on the set of characters of D, denoted |Irr(D)/E|.
# Since D is an elementary abelian 2-group, Irr(D) is isomorphic to D as an E-module.
# So, we need to calculate the number of orbits |D/E|.
# We use Burnside's Lemma: k(B) = (1/|E|) * sum(|D^e| for e in E).

# To use Burnside's Lemma, we need to determine the size of the fixed-point sets |D^e|.
# D is a 5-dimensional vector space over F_2. E = C_5 acts on D.
# The irreducible representations of C_5 over F_2 have dimensions 1 (trivial) and 4.
# Since D is 5-dimensional, its decomposition into irreducible E-modules must be D = V_1 + V_4.
# The fixed-point set D^E corresponds to the trivial submodule V_1, which has dimension 1.
# So, |D^E| = 2^1 = 2.
# For any non-identity element e in E, since E is cyclic, the set of points fixed by e
# is the same as the set of points fixed by the whole group E. So |D^e| = |D^E|.

# For the identity element e=1 in E, all elements of D are fixed.
fixed_points_identity = order_D

# For any non-identity element e in E, the number of fixed points is |D^E|.
dim_fixed_subspace = 1
fixed_points_non_identity = 2**dim_fixed_subspace

# There is 1 identity element and |E|-1 = 4 non-identity elements in E.
num_non_identity_elements = order_E - 1

# Apply Burnside's Lemma
sum_of_fixed_points = fixed_points_identity + (num_non_identity_elements * fixed_points_non_identity)
k_B = sum_of_fixed_points // order_E

print("To calculate k(B), we use Burnside's Lemma for the action of E on D.")
print(f"The identity element of E fixes |D| = {fixed_points_identity} elements.")
print(f"Each of the {num_non_identity_elements} non-identity elements of E fixes {fixed_points_non_identity} elements.")
print(f"By Burnside's Lemma, k(B) = (1/{order_E}) * ({fixed_points_identity} + {num_non_identity_elements} * {fixed_points_non_identity}) = {k_B}.")
print("-" * 20)

# Step 4: Calculate the final result k(B) - l(B).
result = k_B - l_B

print("The final calculation is k(B) - l(B):")
print(f"{k_B} - {l_B} = {result}")
