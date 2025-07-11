import numpy as np
from itertools import combinations

def get_binary_vectors(m):
    """Generates all non-zero binary vectors of dimension m."""
    # There are 2**m possible binary vectors.
    # The first one (all zeros) is excluded as it cannot be part of a linearly independent set.
    vectors = []
    for i in range(1, 2**m):
        bin_str = bin(i)[2:].zfill(m)
        vectors.append(np.array([int(c) for c in bin_str]))
    return vectors

def are_linearly_independent(vectors):
    """Checks if a set of vectors is linearly independent over R."""
    if not vectors:
        return True
    # The vectors are linearly independent if the rank of the matrix
    # formed by these vectors as columns is equal to the number of vectors.
    matrix = np.array(vectors).T
    return np.linalg.matrix_rank(matrix) == len(vectors)

def check_glp(vector_set, m):
    """
    Checks if a set of vectors is in general linear position,
    meaning any subset of size m is linearly independent.
    """
    if len(vector_set) < m:
        return True # Condition is trivially satisfied

    for subset in combinations(vector_set, m):
        if not are_linearly_independent(list(subset)):
            return False
    return True

def find_max_n_by_search(m):
    """
    Finds the maximum n for a given m by exhaustive search.
    This is computationally feasible only for small m.
    """
    all_vectors = get_binary_vectors(m)
    max_n_found = 0
    
    # We test for increasing n
    # Start checking from n=m up to the number of available non-zero vectors
    for n in range(m, len(all_vectors) + 2):
        found_for_n = False
        # If n exceeds available vectors, we stop
        if n > len(all_vectors):
            break
        # Check all combinations of n vectors
        for vector_set in combinations(all_vectors, n):
            if check_glp(list(vector_set), m):
                # print(f"Found a set of n={n} vectors in General Linear Position for m={m}.")
                max_n_found = n
                found_for_n = True
                # Break inner loop and check for the next n
                break 
        
        if not found_for_n:
            # If we couldn't find a set for a given n, we won't find one for n+1.
            # So, the maximum found so far is the overall maximum.
            break
            
    return max_n_found

def solve():
    """
    Solves the problem for a specific m and provides verification.
    """
    # Set m (must be >= 3 as per the problem)
    # Note: Computation for m > 3 will be very slow due to combinatorial explosion.
    m = 3 
    print(f"Determining the maximum value of n for m = {m}")
    print("A set of vectors is in general linear position if any subset of size m is linearly independent.")
    print("-" * 40)

    # 1. Demonstrate that n = m + 1 is possible
    print(f"Testing the constructed set for n = m + 1 = {m+1}...")
    constructed_set = []
    # Standard basis vectors e_i
    for i in range(m):
        v = np.zeros(m, dtype=int)
        v[i] = 1
        constructed_set.append(v)
    # Vector of all ones
    constructed_set.append(np.ones(m, dtype=int))

    is_glp = check_glp(constructed_set, m)
    print("The set tested consists of the standard basis vectors and the vector of all ones:")
    for v in constructed_set:
        print(v)
    if is_glp:
        print(f"\nThis set IS in general linear position. This proves n >= {m+1}.")
    else:
        print(f"\nThis set IS NOT in general linear position.")

    print("-" * 40)
    # 2. Systematically search for the maximum n
    print(f"Searching for the maximum possible n for m={m} by checking all combinations...")
    max_n = find_max_n_by_search(m)
    print("-" * 40)

    print(f"Conclusion for m = {m}:")
    print(f"The maximum value of n found by search is {max_n}.")
    print("The theoretical proof shows the maximum value of n is m+1.")
    
    # Outputting the final equation as requested
    print("\nThe final equation is:")
    print(f"n = m + 1")
    print(f"{max_n} = {m} + 1")

# Run the solution
solve()