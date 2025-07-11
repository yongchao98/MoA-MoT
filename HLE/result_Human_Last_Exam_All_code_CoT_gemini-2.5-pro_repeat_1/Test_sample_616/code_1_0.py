import numpy as np

def solve_brockett_min(A, B):
    """
    Calculates the minimum of the asymmetric Brockett cost function
    f(X, Y) = <A, X^T B Y> for X, Y in SO(n).

    Args:
        A (np.ndarray): The first n x n matrix.
        B (np.ndarray): The second n x n matrix.
    """
    n = A.shape[0]
    if A.shape != (n, n) or B.shape != (n, n):
        print("Error: Matrices must be square and of the same size.")
        return

    print(f"Given matrices:\nA =\n{A}\n\nB =\n{B}\n")

    # Step 1: Compute singular values and sort them in descending order
    a = np.linalg.svd(A, compute_uv=False)
    b = np.linalg.svd(B, compute_uv=False)
    
    # The svd function already returns them sorted descendingly
    print(f"Singular values of A (a_i): {a}")
    print(f"Singular values of B (b_i): {b}\n")

    # Step 2: Create the cross-product terms a_i * b_{n-i+1}
    # In python, b's reversed order is b[::-1]
    b_rev = b[::-1]
    cross_products = a * b_rev
    
    print("Terms to be summed (a_i * b_{n-i+1}):")
    eq_str_parts = []
    for i in range(n):
        eq_str_parts.append(f"{a[i]:.2f} * {b_rev[i]:.2f}")
    print(" + ".join(eq_str_parts))
    print(f"= {cross_products}\n")


    # Step 3: Calculate the two potential minimum values
    sum_val = np.sum(cross_products)
    min_term = np.min(cross_products)

    v1 = -sum_val
    v2 = -sum_val + 2 * min_term
    
    print(f"Base sum S = sum(a_i * b_{n-i+1}) = {sum_val:.4f}")
    print(f"Minimum term min(a_j * b_{n-j+1}) = {min_term:.4f}\n")
    print(f"Potential minimum V1 = -S = {v1:.4f}")
    print(f"Potential minimum V2 = -S + 2 * min_term = {v2:.4f}\n")

    # Step 4: Check the condition based on determinants
    det_A = np.linalg.det(A)
    det_B = np.linalg.det(B)
    sign_det_A = np.sign(det_A)
    sign_det_B = np.sign(det_B)
    
    print("Checking condition:")
    print(f"det(A) = {det_A:.4f}, sign(det(A)) = {sign_det_A}")
    print(f"det(B) = {det_B:.4f}, sign(det(B)) = {sign_det_B}")
    print(f"n = {n}, (-1)^n = {(-1)**n}\n")

    final_min = 0
    reason = ""

    # Check if either matrix is singular
    if np.isclose(det_A, 0) or np.isclose(det_B, 0):
        final_min = v1
        reason = "At least one matrix is singular. Minimum is V1."
    else:
        # Check the sign condition
        if sign_det_A * sign_det_B == (-1)**n:
            final_min = v1
            reason = f"sign(det(A))*sign(det(B)) = {sign_det_A * sign_det_B} == (-1)^n. Minimum is V1."
        else:
            final_min = v2
            reason = f"sign(det(A))*sign(det(B)) = {sign_det_A * sign_det_B} != (-1)^n. Minimum is V2."

    # Step 5: Print the final result and the equation
    print(f"Result:\n{reason}\n")
    
    print("Final Calculation:")
    if final_min == v1:
        print(f"Minimum = -({ ' + '.join([f'{cp:.2f}' for cp in cross_products]) })")
        print(f"Minimum = -({sum_val:.4f}) = {v1:.4f}")
    else:
        print(f"Minimum = -({ ' + '.join([f'{cp:.2f}' for cp in cross_products]) }) + 2 * {min_term:.2f}")
        print(f"Minimum = -({sum_val:.4f}) + {2*min_term:.4f} = {v2:.4f}")

# Example usage with 3x3 matrices
A_example = np.array([[2, 1, 0], [1, 5, 1], [0, 1, 1]])
B_example = np.array([[4, 0, 1], [0, 3, 0], [1, 0, 6]])

solve_brockett_min(A_example, B_example)
