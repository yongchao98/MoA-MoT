import numpy as np

# Step 1: Define the orthonormal basis for C^6.
# These are the first 6 vectors.
basis_vectors = np.identity(6, dtype=complex)
num_basis_vectors = len(basis_vectors)

# Step 2: Define the 4x4 Hadamard matrix to construct orthogonal sets of vectors.
H4 = np.array([
    [1, 1, 1, 1],
    [1, -1, 1, -1],
    [1, 1, -1, -1],
    [1, -1, -1, 1]
])

# Step 3: Define the index sets for the additional vector families.
# Each index set has size 4, and any two sets intersect in 2 elements.
index_sets = [
    [0, 1, 2, 3],  # Corresponds to {e1, e2, e3, e4}
    [0, 1, 4, 5],  # Corresponds to {e1, e2, e5, e6}
    [2, 3, 4, 5]   # Corresponds to {e3, e4, e5, e6}
]

# Step 4: Construct the additional vectors.
# For each index set, create a family of 4 vectors.
additional_vectors_families = []
for indices in index_sets:
    family = []
    for i in range(4):
        new_vector = np.zeros(6, dtype=complex)
        # The new vector is a linear combination of 4 basis vectors.
        # The coefficients are from a row of the Hadamard matrix, scaled by 1/2.
        for j in range(4):
            new_vector[indices[j]] = H4[i, j]
        new_vector *= 0.5
        family.append(new_vector)
    additional_vectors_families.append(family)

num_additional_vectors_per_family = len(additional_vectors_families[0])
num_families = len(additional_vectors_families)

# Step 5: Calculate the total number of vectors.
total_vectors = num_basis_vectors + num_families * num_additional_vectors_per_family

# The final equation representing the total count
print(f"The construction consists of {num_basis_vectors} basis vectors and {num_families} families of {num_additional_vectors_per_family} vectors each.")
print(f"Total number of vectors = {num_basis_vectors} + {num_families} * {num_additional_vectors_per_family}")
# Output each number in the final equation as requested.
print(f"{num_basis_vectors} + {num_additional_vectors_per_family} + {num_additional_vectors_per_family} + {num_additional_vectors_per_family} = {total_vectors}")
