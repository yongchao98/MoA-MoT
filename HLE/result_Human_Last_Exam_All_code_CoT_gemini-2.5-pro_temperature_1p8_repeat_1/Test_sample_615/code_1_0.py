import numpy as np
from itertools import combinations

def check_general_linear_position(vectors, m):
    """
    Checks if a set of vectors is in general linear position.
    This means any subset of size m must form a matrix with a non-zero determinant.
    """
    n = len(vectors)
    if n < m:
        return True

    # Iterate through all combinations of m vectors from the set
    for subset_indices in combinations(range(n), m):
        # Form a matrix from the subset of vectors
        matrix = np.array([vectors[i] for i in subset_indices])
        
        # Check if the vectors are linearly independent by computing the determinant
        if np.isclose(np.linalg.det(matrix), 0):
            print(f"FAIL: Subset with indices {subset_indices} is linearly dependent.")
            print("The vectors are:")
            for i in subset_indices:
                print(f"  {vectors[i]}")
            return False
    
    print("SUCCESS: All tested subsets of size m are linearly independent.")
    return True

# The problem states m >= 3. Let's use m=4 for our demonstration.
m = 4

print("The problem is to find the maximum number n of vectors in {0, 1}^m")
print("such that any m of them are linearly independent.")
print(f"\n--- Testing for n = m + 1 using m = {m} ---\n")

# We construct a set of n = m + 1 vectors.
# The set consists of the m standard basis vectors and the all-ones vector.
n_m_plus_1 = m + 1
construction_m_plus_1 = []
# Add m standard basis vectors e_i
for i in range(m):
    v = np.zeros(m, dtype=int)
    v[i] = 1
    construction_m_plus_1.append(v)
# Add the all-ones vector
construction_m_plus_1.append(np.ones(m, dtype=int))

print(f"For m = {m}, let's test a construction with n = m + 1 = {n_m_plus_1} vectors.")
check_general_linear_position(construction_m_plus_1, m)


print(f"\n--- Testing for n = m + 2 using m = {m} ---\n")

# Now let's try to add one more vector to this set.
# A simple choice for a new binary vector is the sum of the first two basis vectors, e_1 + e_2.
construction_m_plus_2 = list(construction_m_plus_1) # Make a copy
new_vector = np.zeros(m, dtype=int)
new_vector[0] = 1
new_vector[1] = 1
construction_m_plus_2.append(new_vector)
n_m_plus_2 = len(construction_m_plus_2)

print(f"For m = {m}, we test a set with n = m + 2 = {n_m_plus_2} vectors.")
print(f"(The previous {n_m_plus_1} vectors plus the vector {new_vector})")
check_general_linear_position(construction_m_plus_2, m)
# The dependency 1*(e_1+e_2) - 1*(1) + 1*(e_3) + ... + 1*(e_m) = 0 shows these vectors are dependent.
# Our check will find this subset: {e_1+e_2, 1, e_3, ..., e_m}.

print("\n--- Conclusion ---\n")
print("The demonstration shows that n=m+1 is possible, but our attempt")
print("to build a set for n=m+2 failed as it created a linearly dependent subset of size m.")
print("This aligns with the known mathematical result that the maximum is m+1.")
print("\nThe final equation for the maximum value of n is:")
# Output the final relationship as requested
print("n = m + 1")
