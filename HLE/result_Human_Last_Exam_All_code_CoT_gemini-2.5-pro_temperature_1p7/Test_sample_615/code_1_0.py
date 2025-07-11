import numpy as np
from itertools import combinations

def solve_for_m(m):
    """
    Finds the maximum value of n for a given m by exhaustive search.
    """
    print(f"Starting exhaustive search for m = {m}...")

    # Generate all non-zero binary vectors of dimension m.
    # The zero vector is excluded as it would make any set containing it linearly dependent.
    all_vectors = []
    for i in range(1, 2**m):
        vec = tuple(int(c) for c in bin(i)[2:].zfill(m))
        all_vectors.append(vec)

    total_vectors = len(all_vectors)

    # Search downwards for n, starting from the maximum possible number of vectors.
    for n in range(total_vectors, m, -1):
        # Iterate through all combinations of n vectors.
        for v_set in combinations(all_vectors, n):
            is_glp = True  # Assume General Linear Position for this set.
            
            # Check every subset of size m within the current v_set.
            for m_subset in combinations(v_set, m):
                # The vectors are columns in the matrix for independence check.
                matrix = np.array(m_subset).T
                
                # Determinant must be non-zero for linear independence.
                # Since the matrix contains only integers, the determinant is an integer.
                # Using round() for safety with floating point results.
                if round(np.linalg.det(matrix)) == 0:
                    is_glp = False
                    break  # This subset is dependent, so v_set is invalid.
            
            if is_glp:
                # A valid set of size n is found. Since we search downwards, this is the maximum.
                print(f"Search complete. Found maximum n = {n} for m = {m}.")
                return n
    
    # This part should not be reached if m>=1, as n=m is always possible (the basis vectors).
    return m

# For the given problem m >= 3. Let's run for m=3.
m_test = 3
max_n = solve_for_m(m_test)

print(f"\nFor the specific case m = {m_test}, the maximum value of n is {max_n}.")
print("This supports the general theoretical result that n = m + 1.")

# As requested, output the final equation with the numbers from our test case.
print("\nThe final equation is n = m + 1")
print("Using the numbers from our verification:")
# outputting each number in the final equation
print(f"{max_n} = {m_test} + 1")
