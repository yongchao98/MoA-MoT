import numpy as np
from itertools import combinations

def check_general_linear_position(m):
    """
    Constructs and verifies a set of m+1 vectors in {0,1}^m
    that are in general linear position.
    
    Args:
        m (int): The dimension of the vector space, must be >= 3.
    """
    if not isinstance(m, int) or m < 3:
        print("Please provide an integer m >= 3.")
        return

    print(f"Testing for m = {m}")
    
    # 1. Construct the set of n = m + 1 vectors
    vectors = []
    # Add the m standard basis vectors
    for i in range(m):
        v = np.zeros(m, dtype=int)
        v[i] = 1
        vectors.append(v)
    
    # Add the vector of all ones
    vectors.append(np.ones(m, dtype=int))
    
    n = len(vectors)
    print(f"Constructed n = {n} vectors in {{0, 1}}^{m}:")
    for i, v in enumerate(vectors):
        print(f"v_{i+1}: {v}")
        
    # 2. Verify that any subset of m vectors is linearly independent
    all_subsets_independent = True
    
    # Get all combinations of m vectors from the set of n
    indices = range(n)
    for subset_indices in combinations(indices, m):
        # Form a matrix from the subset of vectors
        matrix = np.array([vectors[i] for i in subset_indices]).T
        
        # Calculate the determinant. If it's non-zero, the vectors are independent.
        det = np.linalg.det(matrix)
        
        if np.isclose(det, 0):
            all_subsets_independent = False
            print(f"\nFound a linearly DEPENDENT subset of {m} vectors (indices {tuple(i+1 for i in subset_indices)}):")
            print(matrix)
            print(f"Determinant is {det}")
            break # No need to check further

    print("\n--- Verification Result ---")
    if all_subsets_independent:
        print(f"Success! All C({n}, {m}) = {len(list(combinations(indices, m)))} subsets of size {m} are linearly independent.")
    else:
        print("Failure! Not all subsets of size {m} are linearly independent.")

    # 3. State the final conclusion
    print("\n--- Conclusion ---")
    print("The maximum value of n is m + 1.")
    # The prompt asks to output each number in the final equation.
    print(f"For this case, the maximum n is {m} + 1 = {m+1}")


# You can change this value to test for different m >= 3
m_value = 4
check_general_linear_position(m_value)