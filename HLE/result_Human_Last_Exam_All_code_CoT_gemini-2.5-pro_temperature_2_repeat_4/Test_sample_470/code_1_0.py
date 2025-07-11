import math

# Step 1: Define the parameters from the problem statement.
# D = (C_2)^5, so the order of the defect group D is 2^5.
dim_D = 5
order_D = 2**dim_D

# The inertial quotient E has order 5.
order_E = 5

# Step 2: Calculate l(B), the number of irreducible Brauer characters.
# l(B) is the number of 2'-conjugacy classes, which equals the order of the
# inertial quotient E, as E is abelian and has an order prime to 2.
l_B = order_E

# Step 3: Calculate k(B), the number of irreducible ordinary characters.
# k(B) is the number of conjugacy classes of the semidirect product D x E.
# The formula is k(B) = (1/|E|) * sum(|C_D(g)| for g in E).

# For the identity element e in E, the centralizer C_D(e) is D itself.
size_centralizer_identity = order_D

# For any non-identity element g in E, the action of g on D (a 5-dim vector space
# over F_2) has a fixed-point subspace of dimension 1.
# So, the size of the centralizer is 2^1 = 2.
size_centralizer_non_identity = 2

# There are |E|-1 non-identity elements in E.
num_non_identity_elements = order_E - 1

# Calculate the sum of the sizes of the centralizers.
sum_of_centralizers = size_centralizer_identity + num_non_identity_elements * size_centralizer_non_identity

# Calculate k(B).
# Note: The result of the division must be an integer for the theory to be correct.
k_B = int(sum_of_centralizers / order_E)


# Step 4: Compute and print the final result k(B) - l(B).
difference = k_B - l_B

print(f"The number of irreducible characters, k(B), is {k_B}.")
print(f"The number of Brauer characters, l(B), is {l_B}.")
print(f"The value of k(B) - l(B) is {k_B} - {l_B} = {difference}.")
