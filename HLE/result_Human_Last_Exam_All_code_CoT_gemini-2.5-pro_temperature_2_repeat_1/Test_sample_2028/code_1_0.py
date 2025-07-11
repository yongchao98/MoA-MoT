import numpy as np
import math

def solve_vest_b(S, v, T_matrices, k):
    """
    Solves the VEST problem for the specific case in question (b).
    
    Args:
        S (np.ndarray): The S matrix.
        v (np.ndarray): The v vector.
        T_matrices (list of np.ndarray): The list of T_i matrices.
        k (int): The parameter k.
        
    Returns:
        np.ndarray: The resulting vector.
    """
    n = v.shape[0]
    m = len(T_matrices)
    
    # Classify T_matrices
    # C[p] will store the list of indices i where T_i = E_{p+1, p+1}
    C = [[] for _ in range(n)]
    for i in range(m):
        t = T_matrices[i]
        # Check if t is a matrix with a single 1 on the diagonal
        if np.count_nonzero(t) == 1 and np.sum(t) == 1:
            diag = np.diag(t)
            if np.count_nonzero(diag) == 1:
                p = np.argmax(diag)
                C[p].append(i)

    # c[p] stores |C_p|
    c_counts = [len(Cp) for Cp in C]
    
    print("Plan: Compute U = k! * S * (sum_{p=0 to n-1} C(|C_p|, k) * E_{p+1,p+1}) * v")
    print(f"k = {k}, k! = {math.factorial(k)}")
    print(f"v = {v.flatten()}")
    print("S = \n", S)
    
    total_matrix_sum = np.zeros((n, n))
    
    for p in range(n):
        cp = c_counts[p]
        print(f"\nProcessing p={p} (for basis vector e_{p+1}):")
        print(f"Indices of T_i = E_{{{p+1},{p+1}}}: {C[p]}")
        print(f"Count c_{p+1} = {cp}")
        
        if cp < k:
            binom_val = 0
            print(f"c_{p+1} < k, so C({cp}, {k}) = 0")
        else:
            binom_val = math.comb(cp, k)
            print(f"C({cp}, {k}) = {binom_val}")

        E_p_p = np.zeros((n, n))
        E_p_p[p, p] = 1
        
        term_matrix = binom_val * E_p_p
        total_matrix_sum += term_matrix
        
        if binom_val > 0:
            print(f"Matrix term for p={p}: C({cp}, {k}) * E_{{{p+1},{p+1}}} = {binom_val} * \n{E_p_p}")

    print("\nCalculating final sum matrix part:")
    print(f"Sum Matrix = sum_{{p=0..n-1}} C(c_p,k) * E_p,p = \n{total_matrix_sum}")

    final_matrix = math.factorial(k) * total_matrix_sum
    print(f"Final Matrix = k! * Sum Matrix = {math.factorial(k)} * \n{total_matrix_sum} = \n{final_matrix}")

    result = S @ final_matrix @ v
    
    # Printing the final equation with numbers
    print("\nFinal Equation Breakdown:")
    # We can write S * (k! * sum_matrix) * v
    term1 = S
    term2 = final_matrix
    term3 = v
    print(f"U = S * (k! * sum_matrix) * v")
    print(f"U = \n{term1}\n * \n{term2}\n * \n{term3}")

    print("\nFinal Result:")
    print(f"U = {result.flatten()}")
    return result

# Example from the thinking process
n_val = 4
m_val = 5
k_val = 2
S_mat = np.identity(n_val)
v_vec = np.array([[1], [2], [3], [4]])
T_mats = [
    np.diag([1, 0, 0, 0]),  # p=0
    np.diag([0, 1, 0, 0]),  # p=1
    np.diag([1, 0, 0, 0]),  # p=0
    np.diag([0, 0, 1, 0]),  # p=2
    np.diag([0, 1, 0, 0])   # p=1
]

solve_vest_b(S_mat, v_vec, T_mats, k_val)