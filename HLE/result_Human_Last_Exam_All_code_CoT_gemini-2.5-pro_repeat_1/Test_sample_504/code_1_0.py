# Plan:
# 1. Start with the 6 standard basis vectors in R^6.
# 2. Add vectors based on three special 4-element subsets (supports) of the basis indices.
#    The supports are chosen such that their pairwise intersections are of even size.
#    Let the indices be {1, 2, 3, 4, 5, 6}.
#    The chosen supports are:
#    S1 = {1, 2, 3, 4}
#    S2 = {1, 2, 5, 6}
#    S3 = {3, 4, 5, 6}
#    |S1 intersect S2| = 2, |S1 intersect S3| = 2, |S2 intersect S3| = 2.
# 3. For each support, we can construct 4 vectors using sign patterns from a Hadamard matrix of order 4.
#    These 4 vectors are pairwise orthogonal.
# 4. The total number of vectors is the sum of the basis vectors and the vectors from each support.

num_basis_vectors = 6
num_supports = 3
num_vectors_per_support = 4

total_vectors = num_basis_vectors + num_supports * num_vectors_per_support

print(f"The number of basis vectors is {num_basis_vectors}.")
print(f"We construct vectors on {num_supports} different supports.")
print(f"For each support, we can create {num_vectors_per_support} vectors.")
print(f"The total number of vectors is the sum: {num_basis_vectors} + {num_supports} * {num_vectors_per_support} = {total_vectors}")

# The final result is an integer.
print("\nThe largest number of such vectors is:")
print(total_vectors)