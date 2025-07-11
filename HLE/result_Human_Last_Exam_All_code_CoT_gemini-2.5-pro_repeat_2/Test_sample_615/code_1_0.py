import numpy as np
from itertools import combinations

def check_glp(vectors, m, verbose=True):
    """
    Checks if a set of vectors is in general linear position.
    This means any m vectors from the set must be linearly independent.
    """
    n = len(vectors)
    if n < m:
        if verbose:
            print("Number of vectors is less than m, condition is vacuously true.")
        return True
    
    # Check all combinations of m vectors
    for subset_indices in combinations(range(n), m):
        # Form the matrix from the subset of vectors
        subset_vectors = [vectors[i] for i in subset_indices]
        mat = np.array(subset_vectors)
        
        # Calculate the determinant
        det = np.linalg.det(mat)
        
        if np.isclose(det, 0):
            if verbose:
                print(f"Failed: Subset of vectors with indices {subset_indices} is linearly dependent.")
                print("Vectors:")
                for v in subset_vectors:
                    print(v)
                print(f"Determinant is {det}")
            return False
            
    if verbose:
        print("Success! All subsets of size m are linearly independent.")
    return True

def solve_for_m(m):
    """
    Constructs a set of m+1 binary vectors and verifies they are in GLP.
    """
    print(f"--- Checking for m = {m} ---")
    n = m + 1
    print(f"Constructing a set of n = {n} vectors in {{0, 1}}^{m}.")

    # 1. Create the standard basis vectors
    vectors = [list(np.eye(1, m, k=i)[0]) for i in range(m)]
    
    # 2. Add the vector of all ones
    ones_vector = [1] * m
    vectors.append(ones_vector)
    
    print("The constructed set of vectors is:")
    for v in vectors:
        print(v)
    
    print("\nVerifying if the set is in general linear position...")
    is_in_glp = check_glp(vectors, m)
    
    if is_in_glp:
        print(f"\nA set of n = m + 1 = {n} vectors in {{0, 1}}^{m} in general linear position exists.")
    else:
        print(f"\nThe constructed set of n = {n} is NOT in general linear position.")

# You can change this value to test for different m
m_value = 4
solve_for_m(m_value)

# Based on the analysis, n = m + 1
final_n_value = m_value + 1
# This part is just to satisfy the output format instruction
print(f"\nFor m = {m_value}, the equation for the maximum value n is simply n = m + 1.")
print(f"{final_n_value} = {m_value} + 1")
