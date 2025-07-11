import math

# Step 1: Define the parameters from the problem statement.
# The defect group D = (C_2)^5.
d_order = 2**5
# The inertial quotient E has order 5.
e_order = 5

# Step 2: Calculate l(B), the number of irreducible Brauer characters.
# For a p-block B with inertial quotient E, if p does not divide |E|,
# l(B) is the number of conjugacy classes of E.
# Here p=2 and |E|=5. E is the cyclic group C_5.
# The number of conjugacy classes of an abelian group is its order.
l_B = e_order

# Step 3: Calculate k(B), the number of irreducible ordinary characters.
# For a block with an abelian defect group D, k(B) is the number of
# orbits of the inertial quotient E acting on D.

# We determine the structure of the action of E=C_5 on D=(C_2)^5.
# This is a 5-dim representation of C_5 over F_2. It decomposes into
# a direct sum of a 1-dim trivial and a 4-dim irreducible representation.
# The fixed-point space D^E has dimension 1.
dim_DE = 1
d_E_order = 2**dim_DE

# The number of orbits is the sum of orbits of fixed points and non-fixed points.
# Number of fixed points (orbits of size 1) = |D^E|
# Number of orbits of size 5 = (|D| - |D^E|) / |E|
num_orbits_size_1 = d_E_order
num_orbits_size_5 = (d_order - d_E_order) // e_order
k_B = num_orbits_size_1 + num_orbits_size_5

# Step 4: Calculate the final result and print it.
result = k_B - l_B

print(f"The number of irreducible ordinary characters is k(B) = {k_B}.")
print(f"The number of irreducible Brauer characters is l(B) = {l_B}.")
print(f"The value of k(B) - l(B) is calculated as follows:")
print(f"k(B) - l(B) = {k_B} - {l_B} = {result}")
