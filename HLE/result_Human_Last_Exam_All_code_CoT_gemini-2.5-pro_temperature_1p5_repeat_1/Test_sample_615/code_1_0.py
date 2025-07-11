import numpy as np
from itertools import combinations

def solve_vector_problem():
    """
    This function demonstrates the solution to the general linear position problem
    for a specific case, m=4.
    """
    m = 4
    # According to our derivation, the maximum number of vectors is n = m + 1.
    n = m + 1
    
    print(f"For m = {m}, the maximum value of n is m + 1.")
    print(f"Final Equation: n = m + 1 = {m} + 1 = {n}")
    print("-" * 30)

    # Part 1: Construct the set of n = m+1 vectors and verify the condition.
    print(f"Part 1: Verifying for n = {n}\n")
    
    # The set V consists of the m standard basis vectors and the all-ones vector.
    V = []
    # Add m standard basis vectors
    for i in range(m):
        v = np.zeros(m, dtype=int)
        v[i] = 1
        V.append(v)
    # Add the all-ones vector
    V.append(np.ones(m, dtype=int))

    print("Constructed set of vectors V:")
    for v in V:
        print(v)
    print("\nChecking all subsets of V of size m:")

    is_valid = True
    # Get all combinations of m vectors from V
    for subset_indices in combinations(range(n), m):
        subset_vectors = [V[i] for i in subset_indices]
        
        # Form a matrix from the subset of vectors. The vectors are rows,
        # so we transpose to make them columns for the determinant calculation.
        matrix = np.array(subset_vectors).T
        det = np.linalg.det(matrix)
        
        print(f"Subset: {subset_indices}, Determinant: {det:.1f}")
        if np.isclose(det, 0):
            is_valid = False
            print("  --> This subset is linearly DEPENDENT. Fails.")
        else:
            print("  --> This subset is linearly INDEPENDENT. OK.")
            
    if is_valid:
        print("\nConclusion: The set of n=m+1 vectors is in general linear position.")
    else:
        print("\nConclusion: The set of n=m+1 vectors is NOT in general linear position.")

    print("-" * 30)

    # Part 2: Show that n=m+2 is not possible.
    print(f"Part 2: Showing that n = m + 2 = {m+2} is not possible.\n")
    
    # Let's add one more vector v_new to our set V.
    # v_new must have between 2 and m-1 ones. Let's pick (1,1,0,0).
    v_new = np.array([1, 1, 0, 0])
    
    V_extended = V + [v_new]
    print("Trying to add a new vector:", v_new)

    # We need to find one subset of size m from V_extended that is dependent.
    # Based on our proof, the dependency is v_new = e_0 + e_1.
    # So the subset {v_new, e_0, e_1} is dependent. To make a set of size m=4,
    # we can add any other vector, for example e_2.
    dependent_subset_vectors = [v_new, V[0], V[1], V[2]]
    
    print("\nChecking a specific subset of the extended set: {v_new, e_0, e_1, e_2}")
    for v in dependent_subset_vectors:
        print(v)

    matrix = np.array(dependent_subset_vectors).T
    det = np.linalg.det(matrix)
    print(f"\nDeterminant of this subset: {det:.1f}")
    if np.isclose(det, 0):
        print("  --> This subset is linearly DEPENDENT.")
        print("\nConclusion: Adding a 6th vector creates a dependent subset of size 4.")
        print("Therefore, n cannot be m+2.")
    else:
        print("  --> This subset is linearly INDEPENDENT.")

solve_vector_problem()