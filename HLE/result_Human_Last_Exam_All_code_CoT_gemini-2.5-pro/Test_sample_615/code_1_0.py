import numpy as np
from itertools import combinations

def check_general_linear_position(vectors, m):
    """
    Checks if a set of vectors is in general linear position.
    This means any subset of size m must be linearly independent.
    """
    if len(vectors) < m:
        return True  # Condition is vacuously true
    
    # Generate all subsets of size m
    for subset in combinations(vectors, m):
        # Form a matrix from the subset of vectors (as columns)
        matrix = np.array(subset).T
        # Check if the matrix is singular (determinant is close to zero)
        if np.isclose(np.linalg.det(matrix), 0):
            print(f"  - Subset {list(subset)} is linearly DEPENDENT (det=0).")
            return False
            
    print(f"  + All {len(list(combinations(vectors, m)))} subsets of size {m} are linearly independent.")
    return True

def main():
    m = 3
    print(f"--- Testing for m = {m} ---")

    # 1. Construct the set V for n = m + 1
    print(f"\nPart 1: Testing construction for n = m+1 = {m+1}")
    e_vectors = [tuple(1 if i == j else 0 for j in range(m)) for i in range(m)]
    one_vector = tuple(1 for _ in range(m))
    V = e_vectors + [one_vector]
    
    print("Constructed set V:")
    for v in V:
        print(f"  {v}")
        
    print("\nChecking if V is in general linear position:")
    is_glp = check_general_linear_position(V, m)
    if is_glp:
        print("Conclusion: A set of n=m+1 vectors in general linear position exists.")
    else:
        print("Conclusion: The constructed set is NOT in general linear position.")
        
    # 2. Try to add another vector for n = m + 2
    print(f"\nPart 2: Testing if n can be m+2 = {m+2}")
    
    # Generate all possible binary vectors in {0,1}^m
    all_binary_vectors = [tuple(int(c) for c in f'{i:0{m}b}') for i in range(2**m)]
    
    # Filter out vectors already in V and the zero vector
    candidate_vectors = [v for v in all_binary_vectors if v not in V and sum(v) > 0]
    
    print(f"\nTrying to add one of the remaining {len(candidate_vectors)} vectors to V:")
    
    found_glp_for_m_plus_2 = False
    for v_new in candidate_vectors:
        V_prime = V + [v_new]
        print(f"\nTesting with new vector v = {v_new}")
        
        # Check if this new set of size m+2 is in GLP
        if check_general_linear_position(V_prime, m):
            print(f"Found a set of {m+2} vectors in GLP! This is unexpected.")
            found_glp_for_m_plus_2 = True
            break
        else:
            print(f"-> Adding {v_new} fails the condition.")

    if not found_glp_for_m_plus_2:
        print("\nConclusion: As proven, no vector can be added. Max n is m+1.")

if __name__ == '__main__':
    main()