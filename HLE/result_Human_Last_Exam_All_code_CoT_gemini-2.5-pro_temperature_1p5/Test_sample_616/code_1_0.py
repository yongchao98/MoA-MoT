import numpy as np

def calculate_min_brockett_cost(A, B):
    """
    Calculates the minimum of the asymmetric Brockett cost function
    f(X, Y) = <A, X^T B Y> for X, Y in SO(n).

    Args:
        A (np.ndarray): An n x n matrix.
        B (np.ndarray): An n x n matrix.
    
    Returns:
        float: The minimum value of the cost function.
    """
    # Get the dimension of the matrices
    n = A.shape[0]

    # Check if matrices are square and have the same dimension
    if A.shape != (n, n) or B.shape != (n, n):
        raise ValueError("Matrices A and B must be square and of the same size.")

    # Step 1: Calculate singular values (a_i and b_i)
    # np.linalg.svd returns singular values sorted in descending order.
    _, a, _ = np.linalg.svd(A)
    _, b, _ = np.linalg.svd(B)

    # Step 2: Calculate determinants and the sign factor s(det(A)det(B))
    det_A = np.linalg.det(A)
    det_B = np.linalg.det(B)
    
    # np.sign returns 0 for an input of 0.
    s = np.sign(det_A * det_B)

    # Step 3: Apply the formula for the minimum value
    # The formula is: - (a_1*b_1 + ... + a_{n-1}*b_{n-1}) - s * a_n*b_n
    # Note: Python uses 0-based indexing, so a_k is a[k-1].
    sum_part = np.sum(a[:-1] * b[:-1])
    last_term = s * a[n-1] * b[n-1]
    min_value = -sum_part - last_term

    # Step 4: Print the detailed results
    print("Given matrices:")
    print("A =\n", A)
    print("\nB =\n", B)

    print("\n--- Calculations ---")
    print("Singular values of A (a_i):")
    for i in range(n):
        print(f"a_{i+1} = {a[i]:.4f}")

    print("\nSingular values of B (b_i):")
    for i in range(n):
        print(f"b_{i+1} = {b[i]:.4f}")

    print(f"\ndet(A) = {det_A:.4f}")
    print(f"det(B) = {det_B:.4f}")
    print(f"sign(det(A)*det(B)) = {int(s)}")

    print("\n--- Minimum Value Calculation ---")
    print("The formula for the minimum value is:")
    print("min = - (a_1*b_1 + ... + a_{n-1}*b_{n-1}) - s * a_n*b_n")

    # Print the formula with numbers plugged in
    equation_str = "min = - ("
    for i in range(n - 1):
        equation_str += f"{a[i]:.4f}*{b[i]:.4f}"
        if i < n - 2:
            equation_str += " + "
    equation_str += f") - {int(s)} * {a[n-1]:.4f}*{b[n-1]:.4f}"
    print("\nPlugging in the values:")
    print(equation_str)

    print(f"min = -({sum_part:.4f}) - ({last_term:.4f})")
    print(f"min = {-sum_part:.4f} - {last_term:.4f}")

    print(f"\nThe minimum of the asymmetric Brockett cost function is: {min_value:.4f}")
    return min_value

if __name__ == '__main__':
    # Define two example n x n matrices A and B
    # You can replace these with your own matrices.
    A = np.array([[1, 2, 3], 
                  [4, 9, 6], 
                  [7, 8, 9]], dtype=float)
    
    B = np.array([[9, 8, 1], 
                  [6, 5, 4], 
                  [3, 2, 1]], dtype=float)

    calculate_min_brockett_cost(A, B)
