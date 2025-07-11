# The dimension of any submodule is p(n) = K1*d1(n) + K2*d2(n) + K3*d3(n), where
# d1(n), d2(n), d3(n) are a basis of linearly independent polynomials.
# We need to count the number of distinct coefficient tuples (K1, K2, K3) that can be formed.

# The coefficient K2 = k2, where k2 can be in {0, 1, 2, 3}.
num_k2_choices = 4

# The coefficients K1 and K3 are given by:
# K1 = k1 + k4
# K3 = k3 + k4
# with k1 in {0, 1, 2}, k3 in {0, 1}, and k4 in {0, 1}.

# We will use a set to store the unique pairs (K1, K3).
k1_k3_pairs = set()

# Possible values for k1, k3, k4
k1_range = [0, 1, 2]
k3_range = [0, 1]
k4_range = [0, 1]

# Loop through all combinations to find unique pairs (K1, K3)
for k1 in k1_range:
    for k3 in k3_range:
        for k4 in k4_range:
            K1 = k1 + k4
            K3 = k3 + k4
            k1_k3_pairs.add((K1, K3))

# The number of unique pairs for (K1, K3)
num_k1_k3_pairs = len(k1_k3_pairs)

# The total number of distinct polynomials is the product of the number of choices for K2
# and the number of choices for the pair (K1, K3).
total_polynomials = num_k1_k3_pairs * num_k2_choices

# Output the final calculation
print(f"The number of ways to form the coefficients (K1, K3) is {num_k1_k3_pairs}.")
print(f"The number of ways to form the coefficient K2 is {num_k2_choices}.")
print("The total number of distinct polynomials is the product of these numbers:")
print(f"{num_k1_k3_pairs} * {num_k2_choices} = {total_polynomials}")