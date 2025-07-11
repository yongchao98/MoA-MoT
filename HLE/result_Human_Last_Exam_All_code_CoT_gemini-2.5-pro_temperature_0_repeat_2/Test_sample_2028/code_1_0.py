import math
import numpy as np

def solve_vest_b(n, m, k, v, S, T_indices):
    """
    Solves the VEST problem for the restricted case in question (b).

    Args:
        n (int): Dimension of the vector space.
        m (int): Number of T matrices.
        k (int): The parameter k.
        v (np.array): The vector v.
        S (np.array): The matrix S.
        T_indices (list): A list of length m, where T_indices[i] is the diagonal
                          index (1-based) of the non-zero entry in T_i. 0 if T_i is zero.
    """
    print("--- VEST Problem (b) Calculation ---")
    print(f"n={n}, m={m}, k={k}\n")
    print(f"v = {v}")
    print(f"S = \n{S}\n")

    # Group T matrices by their non-zero index
    M = {j: [] for j in range(1, n + 1)}
    for i, j in enumerate(T_indices):
        if j > 0:
            M[j].append(i + 1)

    print("Step 1: Group T_i matrices.")
    for j in range(1, n + 1):
        print(f"M_{j} (indices of T_i=E_{j,j}): {M[j]}, size = {len(M[j])}")
    print("-" * 20)

    # Pre-calculate v^T * S
    vT_S = v.T @ S
    print(f"Step 2: Calculate v^T * S = {vT_S}")
    print("-" * 20)

    total_sum = 0
    
    print("Step 3: Calculate sum over j=1 to n: C(|M_j|, k) * v_j * (v^T * S)_j")

    # Define a combination function that handles n < r case
    def combinations(n, k):
        if k < 0 or k > n:
            return 0
        return math.comb(n, k)

    for j in range(1, n + 1):
        size_Mj = len(M[j])
        comb = combinations(size_Mj, k)
        v_j = v[j - 1]
        vTS_j = vT_S[j - 1]
        
        term = comb * v_j * vTS_j
        total_sum += term
        
        print(f"Term for j={j}: C({size_Mj}, {k}) * {v_j} * {vTS_j} = {comb} * {v_j} * {vTS_j} = {term}")

    print("-" * 20)
    print(f"Final Result: The total sum is {total_sum}")
    print("--- End of Calculation ---")

if __name__ == '__main__':
    # Example instance based on the thought process
    n_val = 3
    m_val = 5
    k_val = 2
    
    # v = [1, 2, 3]^T
    v_val = np.array([1, 2, 3])
    
    # S = [[1, 0, 1], [0, 1, 0], [1, 0, 1]]
    S_val = np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]])
    
    # T_1=E_11, T_2=E_22, T_3=E_11, T_4=E_33, T_5=E_11
    # We represent T_i by the index of the '1' on the diagonal.
    T_indices_val = [1, 2, 1, 3, 1]
    
    solve_vest_b(n_val, m_val, k_val, v_val, S_val, T_indices_val)