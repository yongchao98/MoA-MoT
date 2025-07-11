import numpy as np

def solve_vest_decision_c(v, S, T_list, k, theta):
    """
    Solves the decision version of VEST for the case in question (c).

    The problem is to determine if value >= theta, where
    value = v^T * S * (sum(T_i))^k * v

    Args:
        v (np.ndarray): The initial n-dimensional vector.
        S (np.ndarray): An n x n matrix.
        T_list (list of np.ndarray): A list of m n x n matrices.
        k (int): The parameter k.
        theta (float): The threshold for the decision problem.

    Returns:
        bool: True if the computed value is >= theta, False otherwise.
    """
    n = v.shape[0]
    m = len(T_list)

    if m == 0:
        T_sum = np.zeros((n, n))
    else:
        # Step 1: Compute T_sum = sum(T_i for i in 1..m)
        T_sum = np.sum(T_list, axis=0)

    # Step 2: Compute T_sum_k_v = (T_sum^k) * v
    # This is done by k matrix-vector multiplications to be efficient.
    # T_sum_k = np.linalg.matrix_power(T_sum, k) # Less efficient for large k
    T_sum_k_v = v
    for _ in range(k):
        T_sum_k_v = T_sum @ T_sum_k_v

    # Step 3: Compute the final value = v^T * S * T_sum_k_v
    final_value = v.T @ S @ T_sum_k_v
    
    # --- Outputting the equation as requested ---
    # NOTE: The full equation can be very long. We show the components.
    print("--- VEST Decision Problem ---")
    print(f"n (dimension) = {n}")
    print(f"m (num matrices) = {m}")
    print(f"k = {k}")
    print(f"theta = {theta}\n")

    print("Input Components:")
    print(f"v:\n{v}\n")
    print(f"S:\n{S}\n")
    # Limiting print output for large number of T matrices
    for i, T in enumerate(T_list[:min(m, 5)]):
        print(f"T_{i+1}:\n{T}\n")
    if m > 5:
        print("... (more T matrices not shown)\n")

    print("Intermediate Calculations:")
    print(f"T_sum = sum(T_i):\n{T_sum}\n")
    # For demonstration, T_sum^k is computed. Not needed for the efficient algorithm.
    try:
        T_sum_k = np.linalg.matrix_power(T_sum, k)
        print(f"T_sum^k:\n{T_sum_k}\n")
    except np.linalg.LinAlgError:
        print("T_sum^k could not be computed (e.g., singular matrix to a negative power).\n")

    print("Final Equation and Result:")
    print("v^T * S * (T_sum)^k * v = ?")
    print(f"Value computed = {final_value}")
    
    result = final_value >= theta
    print(f"Is {final_value} >= {theta}?  {result}\n")

    return result

if __name__ == '__main__':
    # Example usage based on the thought process.
    # n=2, m=2, k=2
    n_val = 2
    v_vec = np.array([1, 1])
    S_mat = np.identity(n_val)
    # T matrices with one non-zero entry per row
    T1 = np.array([[1, 0], [3, 0]], dtype=float)
    T2 = np.array([[0, 5], [4, 0]], dtype=float)
    T_matrices = [T1, T2]
    k_val = 2
    theta_val = 100

    solve_vest_decision_c(v_vec, S_mat, T_matrices, k_val, theta_val)
