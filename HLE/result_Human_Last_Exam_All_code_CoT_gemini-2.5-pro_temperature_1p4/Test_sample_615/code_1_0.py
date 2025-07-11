import numpy as np
import itertools

def find_max_n(m):
    """
    Finds the maximum number (n) of m-dimensional binary vectors in general
    linear position.

    Args:
        m (int): The dimension of the vector space.

    Returns:
        int: The maximum number of vectors, n.
    """
    # Generate all non-zero binary vectors of dimension m
    all_vectors = list(itertools.product([0, 1], repeat=m))
    # The first vector is the zero vector, which cannot be part of any
    # linearly independent set, so we exclude it.
    non_zero_vectors = all_vectors[1:]

    # Iterate downwards from the maximum possible size of a subset
    for n in range(len(non_zero_vectors), m - 1, -1):
        # Check all subsets of size n
        for subset in itertools.combinations(non_zero_vectors, n):
            is_valid_set = True
            # Check all combinations of m vectors within the subset
            for m_vectors in itertools.combinations(subset, m):
                # Form a matrix and check if its determinant is close to zero
                matrix = np.array(m_vectors).T
                if np.isclose(np.linalg.det(matrix), 0):
                    is_valid_set = False
                    break  # This combination of m vectors is not LI
            
            if is_valid_set:
                # Found a valid subset of size n. Since we iterate downwards,
                # this is the maximum possible n.
                return n
    return 0

# Set the dimension m. m=4 disambiguates the answer choices.
m = 4
max_n = find_max_n(m)

print(f"For m = {m}, the maximum value of n is {max_n}.")
print("The relationship is n = m + 1.")
print("So, for m = 4, the equation is:")
# Output each number in the final equation
print(f"{max_n} = {m} + 1")
