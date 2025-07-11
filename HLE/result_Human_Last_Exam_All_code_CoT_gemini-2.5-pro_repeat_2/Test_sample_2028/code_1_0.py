import numpy as np

def solve_vest_a(v, T_list, S, k):
    """
    Solves an instance of the VEST problem under the restrictions of question (a).
    
    The problem is to compute y = S * (sum of T_i)^k * v, which holds for
    commuting matrices when summing over all k-sequences of transformations.

    Args:
        v (np.ndarray): The initial rational vector.
        T_list (list of np.ndarray): A list of rational matrices that commute pairwise.
        S (np.ndarray): A rational matrix. For question (a), this is the identity matrix.
        k (int): The parameter k.
    """
    
    # The dimension of the vector space
    n = v.shape[0]
    
    # The number of transformations
    m = len(T_list)

    # Sanity check: verify that S is the identity matrix as per the problem statement
    if not np.allclose(S, np.identity(n)):
        print("Warning: S is not the identity matrix, but proceeding with calculation.")

    # Sanity check: verify that matrices in T_list commute pairwise.
    # This is computationally expensive and not required for the calculation itself,
    # but good for verifying the input satisfies the problem's constraints.
    # for i in range(m):
    #     for j in range(i + 1, m):
    #         if not np.allclose(T_list[i] @ T_list[j], T_list[j] @ T_list[i]):
    #             raise ValueError("Matrices in T_list do not commute pairwise.")

    # Step 1: Compute the sum of all T_i matrices.
    T_sum = np.zeros_like(T_list[0], dtype=float)
    for T in T_list:
        T_sum += T

    # Step 2: Compute T_sum to the power of k using binary exponentiation (matrix exponentiation).
    # This is an efficient way to compute matrix powers.
    T_sum_k = np.linalg.matrix_power(T_sum, k)

    # Step 3: Compute the final vector y = S * T_sum_k * v
    y = S @ T_sum_k @ v
    
    print("The final equation is of the form: y = S * (T_1 + ... + T_m)^k * v")
    print(f"For the given instance, the resulting vector y has {y.shape[0]} components:")
    for i in range(y.shape[0]):
        print(f"y[{i}] = {y[i]}")

if __name__ == '__main__':
    # This is a sample instance for the VEST problem as described in (a).
    # Let n=2 (dimension), m=2 (number of matrices), k=3.
    
    # S is the identity matrix.
    S_matrix = np.identity(2)
    
    # v is a rational vector.
    v_vector = np.array([1, 2])
    
    # T_i are commuting rational matrices.
    # T1 = [[1, 2], [2, 1]], T2 = [[2, 1], [1, 2]]
    # T1*T2 = [[4, 5], [5, 4]]
    # T2*T1 = [[4, 5], [5, 4]]
    # So they commute.
    T1_matrix = np.array([[1, 2], [2, 1]])
    T2_matrix = np.array([[2, 1], [1, 2]])
    T_matrices = [T1_matrix, T2_matrix]
    
    k_param = 3
    
    # Solve the problem for the sample instance.
    solve_vest_a(v_vector, T_matrices, S_matrix, k_param)
    
    print("\n---")
    print("The answers to the theoretical questions are:")
    print("(a) No")
    print("(b) Yes")
    print("(c) W[1]-hard")

