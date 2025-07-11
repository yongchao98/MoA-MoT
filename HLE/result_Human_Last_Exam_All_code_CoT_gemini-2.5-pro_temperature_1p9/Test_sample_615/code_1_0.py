import numpy as np
from itertools import combinations

def check_general_linear_position(vectors, k):
    """
    Checks if any subset of size k from the list of vectors is linearly independent.
    Returns True if the condition holds, False otherwise.
    Also returns the first found dependent set if any.
    """
    if len(vectors) < k:
        return True, None

    vector_indices = range(len(vectors))
    
    for subset_indices in combinations(vector_indices, k):
        # Form the matrix from the subset of vectors
        # Note: numpy wants rows, but we think of them as columns, so we transpose.
        matrix = np.array([vectors[i] for i in subset_indices]).T
        
        # Check for linear independence by calculating the determinant
        if np.linalg.det(matrix) == 0:
            dependent_set = [vectors[i] for i in subset_indices]
            return False, dependent_set
            
    return True, None

def solve():
    """
    Solves the problem and provides a demonstration.
    """
    # Let's choose m=4 for demonstration, as per the problem m >= 3.
    m = 4
    
    # Based on our derivation, the maximum value of n should be m + 1.
    n_max = m + 1
    
    print(f"Problem: Determine the maximum n for m-dimensional binary vectors in general linear position.")
    print(f"Our derived answer is n = m + 1.\n")

    print(f"--- Demonstration for m = {m} ---")
    print(f"Maximum n should be {n_max}.\n")
    
    # 1. Construct a set of n = m + 1 vectors.
    # The set consists of the m standard basis vectors and the all-ones vector.
    vectors_m_plus_1 = []
    # Add standard basis vectors e_i
    for i in range(m):
        v = [0] * m
        v[i] = 1
        vectors_m_plus_1.append(tuple(v))
    # Add the all-ones vector
    vectors_m_plus_1.append(tuple([1] * m))

    print(f"Let's test the set of n = {n_max} vectors:")
    for v in vectors_m_plus_1:
        print(v)
    print("")

    # 2. Verify that any m vectors from this set are linearly independent.
    print(f"Checking if any m = {m} vectors from this set are linearly independent...")
    is_independent, _ = check_general_linear_position(vectors_m_plus_1, m)

    if is_independent:
        print(f"Result: Success! All subsets of size {m} are linearly independent.")
        print(f"This supports the conclusion that n >= {m+1}.\n")
    else:
        print(f"Result: Failure! The construction is wrong (which would be unexpected).")

    # 3. Show that adding one more vector can break the property.
    print(f"Now, let's try to add a ({m+2})-th vector and see if the property breaks.")
    vectors_m_plus_2 = list(vectors_m_plus_1)
    # Add an arbitrary vector, e.g., (1, 1, 0, 0, ...)
    new_vector = tuple([1, 1] + [0] * (m - 2))
    vectors_m_plus_2.append(new_vector)
    
    print(f"The new set of n = {m+2} vectors is:")
    for v in vectors_m_plus_2:
        print(v)
    print("")

    print(f"Checking if any m = {m} vectors from this new set are linearly independent...")
    is_independent, dependent_set = check_general_linear_position(vectors_m_plus_2, m)

    if not is_independent:
        print(f"Result: Success! We found a linearly DEPENDENT subset of size {m}:")
        for v in dependent_set:
            print(v)
        print(f"This supports the conclusion that n cannot be {m+2}.")
    else:
        print("Could not find a dependent set with this choice. The proof is general, though.")

    # Final equation part of the prompt
    print("\nFinal Conclusion:")
    print(f"The maximum value of n is a function of m.")
    print(f"The final equation for n is: n = m + 1")


solve()