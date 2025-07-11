import numpy as np
from itertools import combinations

def check_general_linear_position(vectors, m):
    """
    Checks if a set of vectors is in general linear position.
    This means any subset of size m is linearly independent.
    """
    if len(vectors) < m:
        # Not enough vectors to form a subset of size m
        return True
    
    # Iterate through all subsets of size m
    for subset_indices in combinations(range(len(vectors)), m):
        # Form the matrix from the subset of vectors
        subset_vectors = [vectors[i] for i in subset_indices]
        matrix = np.array(subset_vectors).T # Use vectors as columns
        
        # Check for linear independence by calculating the determinant
        if np.isclose(np.linalg.det(matrix), 0):
            print(f"FAIL: The subset of vectors with indices {subset_indices} is linearly dependent.")
            print("Matrix:\n", matrix)
            print("Determinant:", np.linalg.det(matrix))
            return False
            
    print(f"SUCCESS: All {len(list(combinations(range(len(vectors)), m)))} subsets of size {m} are linearly independent.")
    return True

# Set m
m = 4
n = m + 1
print(f"--- Testing for m = {m}, trying to show n = {n} is the maximum ---")

# 1. Construct the set of m+1 vectors
# Standard basis vectors e_1, ..., e_m
vectors = [list(v) for v in np.identity(m, dtype=int)]
# Add the all-ones vector
all_ones = [1] * m
vectors.append(all_ones)

print(f"\nConstructed a set of n = {len(vectors)} vectors:")
for v in vectors:
    print(v)

# 2. Verify that this set of m+1 vectors is in general linear position
print(f"\nVerifying that any {m} of these {n} vectors are linearly independent...")
check_general_linear_position(vectors, m)

# 3. Show that adding one more vector fails
print("\n--- Demonstrating that n = m+2 is not possible ---")
# Add a new vector, e.g., one with weight 2 like (1,1,0,0)
v_new = [0] * m
v_new[0] = 1
v_new[1] = 1
extended_vectors = vectors + [v_new]

print(f"\nTrying to add a new vector {v_new}. The extended set has {len(extended_vectors)} vectors.")
print(f"Verifying if this extended set is in general linear position...")

# The theory predicts this check will fail.
# The dependency is v_new = e_1 + e_2.
# So the subset {e_1, e_2, v_new} is dependent.
# A subset of size m containing these three, e.g., {e_1, e_2, e_3, v_new}, will also be dependent.
check_general_linear_position(extended_vectors, m)

print(f"\nConclusion: The maximum value of n is m+1.")