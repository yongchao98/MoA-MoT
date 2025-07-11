# The vector space is C^6. We want to find the largest set of vectors
# where the angle between any two is pi/2 or pi/3, and at least one pair is orthogonal.
# This means for any two distinct unit vectors v, w in the set, |(v,w)| is 0 or 1/2.

# Our construction combines the standard basis vectors with three families of "special" vectors.

# 1. The standard basis vectors. There are 6 of them in C^6.
# They are all mutually orthogonal, so this part of the set satisfies the conditions.
num_basis_vectors = 6

# 2. Special vectors. These vectors are supported on specific subsets of indices of size 4.
# The indices {1, ..., 6} are partitioned into 3 pairs: {1,2}, {3,4}, {5,6}.
# The support sets are formed by taking unions of two of these pairs. There are 3 such support sets.
num_supports = 3

# 3. For each support set of size 4, we define vectors by assigning signs to the 4 components.
# A vector has components +-1/2 on its support and 0 elsewhere.
# There are 2^4 = 16 sign combinations.
# To ensure the vectors are not parallel (angle 0), we choose one vector from each of the 8 pairs {v, -v}.
# This gives 8 vectors per support set that meet the angle criteria among themselves.
num_special_vectors_per_support = 2**4 // 2

# The total number of vectors is the sum from these sources.
total_vectors = num_basis_vectors + num_supports * num_special_vectors_per_support

# We print out each number that goes into the final calculation.
print(f"The total number of vectors is the sum of basis vectors and special vectors from each support set.")
print(f"Number of basis vectors: {num_basis_vectors}")
print(f"Number of support sets for special vectors: {num_supports}")
print(f"Number of special vectors per support set: {num_special_vectors_per_support}")
print(f"The final calculation is: {num_basis_vectors} + {num_supports} * {num_special_vectors_per_support} = {total_vectors}")