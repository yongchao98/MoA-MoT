import numpy as np
from itertools import combinations

def check_general_linear_position(m):
    """
    Constructs a set of m+1 binary vectors and verifies they are in
    general linear position for a given m >= 3.
    """
    if m < 3:
        print("This problem is for m >= 3.")
        return

    n = m + 1
    print(f"For m = {m}, we propose a set of n = m + 1 = {n} vectors.")

    # 1. Construct the set of n vectors
    vectors = []
    # Add the m standard basis vectors
    for i in range(m):
        v = np.zeros(m, dtype=int)
        v[i] = 1
        vectors.append(v)

    # Add the vector of all ones
    vectors.append(np.ones(m, dtype=int))

    print("Constructed Vectors:")
    for i, v in enumerate(vectors):
        print(f"  v{i+1}: {v}")

    # 2. Generate all combinations of m vectors from the set
    subsets = list(combinations(range(n), m))
    print(f"\nChecking all {len(subsets)} subsets of size {m} for linear independence:")

    all_independent = True
    for i, subset_indices in enumerate(subsets):
        # 3. For each subset, form a matrix and calculate its determinant
        matrix_list = [vectors[j] for j in subset_indices]
        matrix = np.array(matrix_list).T # Transpose to have vectors as columns

        # Use round to handle potential floating point inaccuracies
        determinant = round(np.linalg.det(matrix))

        subset_v_names = "{" + ", ".join([f"v{j+1}" for j in subset_indices]) + "}"
        print(f"\nSubset {i+1}: {subset_v_names}")
        print("Matrix:")
        print(matrix)
        print(f"Determinant = {determinant}")

        if determinant == 0:
            print("This subset is linearly DEPENDENT. The condition fails.")
            all_independent = False
            break
        else:
            print("This subset is linearly INDEPENDENT.")

    print("\n" + "="*30)
    if all_independent:
        print(f"Conclusion: The construction for n = {n} is valid as all subsets of size {m} are linearly independent.")
    else:
        print(f"Conclusion: The construction for n = {n} is invalid.")
        
    print("\nThe theoretical maximum value for n is m+1.")


# Run the demonstration for a specific value of m, for example m=4
# For m=4, the maximum value of n is 5.
# This contradicts choice D (n = 2^(4-1) = 8) and confirms choice B (n = 4+1 = 5).
m_value = 4
check_general_linear_position(m_value)
