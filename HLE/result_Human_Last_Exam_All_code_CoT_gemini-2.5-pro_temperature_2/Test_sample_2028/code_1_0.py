import numpy as np

def solve_vest(w, v, T_list, k):
    """
    Calculates the VEST problem value using an efficient FPT algorithm.
    The value is w^T * (sum_{s in [m]^k} T_{s_k} ... T_{s_1}) * v.
    The core calculation relies on the identity:
    sum_{s in [m]^k} T_{s_k} ... T_{s_1} = ( (sum_{i=1 to m} T_i^T)^k )^T
    """
    if not T_list:
        print("The list of matrices is empty.")
        return 0

    n = T_list[0].shape[0]
    print(f"--- VEST Problem Instance ---")
    print(f"Parameter k = {k}")
    print(f"Dimension n = {n}")
    print(f"Vector w:\n{w.flatten()}")
    print(f"Vector v:\n{v.flatten()}")
    print("\nTransformation Matrices T_i:")
    for i, T in enumerate(T_list):
        print(f"T_{i+1}:\n{T}")

    # Step 1: Compute A_T = sum of transposed T_i matrices
    # A_T = (T_1 + T_2 + ... + T_m)^T = T_1^T + T_2^T + ... + T_m^T
    A_T = np.zeros_like(T_list[0], dtype=np.float64)
    for T in T_list:
        A_T += T.T
    
    print("\n--- Calculation Steps ---")
    print(f"Step 1: Compute A_T = sum(T_i^T)\n{A_T}")

    # Step 2: Compute (A_T)^k using exponentiation by squaring (FPT in k)
    A_T_k = np.linalg.matrix_power(A_T, k)
    print(f"\nStep 2: Compute (A_T)^k\n{A_T_k}")

    # Step 3: Transpose the result to get the final matrix M
    M = A_T_k.T
    print(f"\nStep 3: Compute M = ((A_T)^k)^T\n{M}")

    # Step 4: Calculate the final result w^T * M * v
    result = w.T @ M @ v
    print(f"\nStep 4: Final Result = w^T * M * v")
    # To show the final equation, we print the components
    print(f"{w.T.flatten()} * \n{M} * \n{v.flatten()}")
    
    return result[0, 0]

if __name__ == '__main__':
    # Create a sample instance for VEST compliant with question (c)
    # (T_i have only one non-zero entry in each row)
    n = 3
    k = 2
    
    # Define vectors w and v
    w = np.array([[1], [2], [3]])
    v = np.array([[4], [5], [6]])

    # Define matrices T_1 and T_2 with one non-zero entry per row
    T1 = np.array([[0, 2, 0], 
                   [1, 0, 0], 
                   [0, 0, 3]], dtype=np.float64)
    
    T2 = np.array([[5, 0, 0], 
                   [0, 0, 4], 
                   [0, 6, 0]], dtype=np.float64)
    
    T_list = [T1, T2]

    # The existence of this FPT algorithm provides the answers.
    # The complexity is dominated by matrix exponentiation, which is O(n^3 log k),
    # and forming the sum matrix, which is O(m * n^2). This is FPT in k.
    #
    # (a) FPT -> not #W[2]-hard (unless FPT=#W[2]). Answer: No.
    # (b) FPT -> not #W[1]-hard (unless FPT=#W[1]). Answer: No.
    # (c) The complexity of the decision problem is determined by the complexity
    #     of computing the value, which is FPT. Answer: FPT.

    final_value = solve_vest(w, v, T_list, k)

    print("\nFinal calculated value:", final_value)
    
    print("\n--- Answers to the questions ---")
    # Directly providing the final answers as requested.
    print("(a) [No]; (b) [No]; (c) [FPT]")


<<< (a) [No]; (b) [No]; (c) [FPT] >>>