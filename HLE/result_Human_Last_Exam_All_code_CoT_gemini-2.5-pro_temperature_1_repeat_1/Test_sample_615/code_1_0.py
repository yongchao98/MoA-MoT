import numpy as np
from itertools import combinations

def check_general_linear_position(vectors, m):
    """
    Checks if a set of vectors in {0, 1}^m is in general linear position.
    This means any subset of m vectors is linearly independent.
    """
    n = len(vectors)
    if n < m:
        # Not enough vectors to form a subset of size m
        print(f"Number of vectors ({n}) is less than the dimension ({m}).")
        return True

    # Iterate through all subsets of size m
    for subset_indices in combinations(range(n), m):
        # Form the matrix from the subset of vectors
        matrix = np.array([vectors[i] for i in subset_indices]).T
        
        # Check for linear independence by computing the determinant
        det = np.linalg.det(matrix)
        
        if np.isclose(det, 0):
            print(f"Subset of vectors at indices {subset_indices} is linearly dependent.")
            print("The vectors are:")
            for i in subset_indices:
                print(vectors[i])
            print(f"Matrix:\n{matrix}")
            print(f"Determinant: {det}")
            return False
            
    print(f"All subsets of size {m} are linearly independent.")
    return True

# Let's test for m=3. The maximum n should be m+1 = 4.
m = 3
n = m + 1

# Our constructed set of vectors for n = m+1
vectors = []
# Add the standard basis vectors e_1, ..., e_m
for i in range(m):
    vec = [0] * m
    vec[i] = 1
    vectors.append(tuple(vec))
# Add the vector of all ones
vectors.append(tuple([1] * m))

print(f"For m = {m}, testing a set of n = {n} vectors:")
for v in vectors:
    print(v)

check_general_linear_position(vectors, m)

# Let's try to add a 5th vector for m=3, n=5
# Any other non-zero vector in {0,1}^3 is a sum of two basis vectors.
# e.g., (1,1,0) = (1,0,0) + (0,1,0)
vectors_n5 = list(vectors)
vectors_n5.append((1,1,0))
print("\n" + "="*20)
print(f"For m = {m}, testing a set of n = {n+1} vectors:")
for v in vectors_n5:
    print(v)
# This will fail because {(1,0,0), (0,1,0), (1,1,0)} is a linearly dependent subset of size 3.
check_general_linear_position(vectors_n5, m)

# Final conclusion based on the mathematical proof
print(f"\nConclusion: For m >= 3, the maximum value of n is m + 1.")
final_n = "m+1"
print(f"The maximum value of n is {final_n}")
