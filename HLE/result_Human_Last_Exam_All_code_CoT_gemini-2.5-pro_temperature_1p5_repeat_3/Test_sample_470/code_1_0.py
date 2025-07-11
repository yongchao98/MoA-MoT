import math

# Step 1: Define the given parameters.
# D = (C_2)^5, so p=2 and the dimension of D as a vector space over F_2 is 5.
dim_D = 5
D_size = 2**dim_D

# The inertial quotient E has order 5.
E_order = 5

# Step 2: Calculate l(B), the number of irreducible Brauer characters.
# For a block with an abelian defect group and p'-inertial quotient E,
# l(B) is the number of irreducible characters of E of p-defect zero.
# Since p=2 and |E|=5 are coprime, all irreducible characters of E have defect zero.
# As E is abelian, the number of its irreducible characters is its order.
l_B = E_order

# Step 3: Calculate k(B), the number of ordinary characters, using Burnside's Lemma.
# k(B) is the number of orbits of E on Irr(D).
# |Irr(D)| = |D| = 32.
# Number of orbits = (1/|E|) * sum_{g in E} |Fixed points of g|.

# For the identity element e in E, all characters are fixed.
num_fixed_by_id = D_size

# For any non-identity element g in E, we determine the number of fixed points.
# The action of g on D=(C_2)^5 corresponds to a linear map of order 5 on a 5D vector space over F_2.
# For a non-trivial action, this space decomposes into a 1D trivial representation
# (the fixed-point space) and a 4D irreducible representation.
# The dimension of the fixed-point space is 1.
# Number of fixed points is 2^1 = 2.
dim_fixed_space = 1
num_fixed_by_gen = 2**dim_fixed_space

# There is 1 identity element and (E_order - 1) non-identity elements.
sum_of_fixed_points = num_fixed_by_id + (E_order - 1) * num_fixed_by_gen

# Apply Burnside's Lemma
k_B = sum_of_fixed_points / E_order

# Ensure k_B is an integer.
k_B = int(k_B)

# Step 4: Compute the final value k(B) - l(B).
result = k_B - l_B

# Output the results clearly.
print(f"The number of irreducible ordinary characters, k(B), is {k_B}.")
print(f"The number of irreducible Brauer characters, l(B), is {l_B}.")
print(f"The calculation for k(B) - l(B) is: {k_B} - {l_B} = {result}")

print(f"<<<{result}>>>")