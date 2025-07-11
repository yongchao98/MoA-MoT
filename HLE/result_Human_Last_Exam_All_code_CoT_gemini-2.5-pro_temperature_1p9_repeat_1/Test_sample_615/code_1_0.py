import numpy as np
from itertools import combinations

def solve_vector_problem():
    """
    This script demonstrates that for a given dimension m, it is possible to find
    n = m + 1 binary vectors in general linear position.
    It does so by constructing a specific set of vectors and verifying the property.

    The problem asks for the maximum n. The theoretical maximum is n = m + 1.
    Proof Outline:
    1. Proof that n <= m + 1: Assume n = m + 2 vectors exist. Pick any m as a basis.
       Express v_{m+1} and v_{m+2} in this basis. This leads to a contradiction
       when you test the linear independence of the set {v_1,...,v_{m-1}, v_{m+1}, v_{m+2}}.
    2. Proof that n = m + 1 is achievable (Constructive Proof):
       The set consisting of the m standard basis vectors {e_1, ..., e_m} plus
       the vector of all ones, j = (1, ..., 1), satisfies the condition. Any
       subset of m vectors from this set is linearly independent.
    
    This script implements the constructive proof for a specific m.
    """
    
    # Let's choose m >= 3 as per the problem. Let's use m=4 for demonstration.
    m = 4
    n = m + 1
    
    print(f"Let m = {m}. We will show that n = m + 1 = {n} is achievable.\n")

    # Construct the set of n = m + 1 vectors
    vectors = []
    # Add the m standard basis vectors
    for i in range(m):
        v = np.zeros(m, dtype=int)
        v[i] = 1
        vectors.append(v)
    
    # Add the vector of all ones
    vectors.append(np.ones(m, dtype=int))

    print(f"We construct the following set of {n} vectors in {{0, 1}}^{m}:")
    for i, v in enumerate(vectors):
        print(f"v_{i+1} = {list(v)}")

    print(f"\nNow, we verify that any subset of {m} vectors is linearly independent.")
    print("-" * 60)

    # Get all combinations of m vectors from the set
    vector_indices = range(n)
    subsets_indices = combinations(vector_indices, m)

    all_independent = True
    for i, subset_idx in enumerate(subsets_indices):
        subset_vectors = [vectors[j] for j in subset_idx]
        
        # Form a matrix from the subset of vectors (as columns)
        matrix = np.array(subset_vectors).T
        
        # Calculate the determinant
        det = np.linalg.det(matrix)
        
        print(f"Subset #{i+1}, consisting of vectors {tuple(j+1 for j in subset_idx)}:")
        # Pretty print the matrix
        matrix_str = str(matrix).replace(' [', '[[').replace(']]', ']]').replace('\n ', ',\n        ')
        
        print(f"  Matrix = {matrix_str}")
        print(f"  det({list(list(v) for v in subset_vectors)}) = {det:.1f}")

        if np.isclose(det, 0):
            print("  --> This set is linearly DEPENDENT. The construction failed.")
            all_independent = False
            break
        else:
            print("  --> The determinant is non-zero, so the vectors are linearly INDEPENDENT.")
        print("-" * 60)

    if all_independent:
        print("\nConclusion:")
        print("All subsets of size m are linearly independent.")
        print(f"This demonstrates that n = m + 1 is achievable.")
        print("Combined with the proof that n cannot exceed m + 1, the maximum value is n = m + 1.")

solve_vector_problem()
<<<B>>>