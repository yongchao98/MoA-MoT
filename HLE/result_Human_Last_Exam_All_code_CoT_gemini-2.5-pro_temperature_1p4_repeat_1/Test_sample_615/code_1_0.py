import numpy as np
import itertools

def find_max_n_for_m(m):
    """
    Finds the maximum n for a given m such that there exist n binary vectors 
    in {0, 1}^m in general linear position.
    
    A set of vectors is in general linear position if any subset of m vectors 
    is linearly independent.
    """
    if m < 3:
        print(f"The problem is stated for m >= 3. Testing for m={m} is not requested.")
        return

    # 1. Generate all non-zero binary vectors of length m.
    all_vectors_tuples = list(itertools.product([0, 1], repeat=m))
    all_vectors_tuples.remove(tuple([0] * m))
    
    # Convert tuples to lists for easier handling with numpy
    all_vectors = [list(v) for v in all_vectors_tuples]
    
    num_total_vectors = len(all_vectors)

    # 2. Search for the largest n downwards.
    max_n_found = 0
    for n in range(num_total_vectors, m - 1, -1):
        # 3. Iterate through all subsets of vectors of size n.
        for vector_subset in itertools.combinations(all_vectors, n):
            is_in_glp = True
            # 4. For each subset, check the general linear position property.
            #    This means every subset of size m must be linearly independent.
            for m_subset in itertools.combinations(vector_subset, m):
                # The vectors in m_subset are the rows of the matrix we want to test.
                # We need to transpose it so vectors are columns.
                matrix = np.array(m_subset).T
                
                # 5. Check for linear independence using the determinant.
                if np.linalg.det(matrix) == 0:
                    is_in_glp = False
                    break  # This m_subset is dependent, so vector_subset is not valid.
            
            if is_in_glp:
                # 6. If a valid set is found, this n is the maximum.
                max_n_found = n
                print(f"For m = {m}, the maximum value of n is {max_n_found}.")
                print("This corresponds to the equation: n = m + 1")
                print(f"Substituting the values: {max_n_found} = {m} + 1")
                return

    if max_n_found == 0:
        print(f"For m = {m}, no set satisfying the condition was found.")

if __name__ == '__main__':
    # Run the analysis for m = 3
    find_max_n_for_m(3)
    print("-" * 20)
    # Run the analysis for m = 4 to confirm the pattern
    find_max_n_for_m(4)