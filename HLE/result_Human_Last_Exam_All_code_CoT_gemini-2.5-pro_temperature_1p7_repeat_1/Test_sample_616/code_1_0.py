import numpy as np

def solve_brockett_minimum():
    """
    Calculates the minimum of the asymmetric Brockett cost function for two matrices A and B.
    The function f(X, Y) = <A, X^T B Y> for X, Y in SO(n).
    """
    # Define two n x n matrices A and B.
    # You can change these matrices to your own.
    A = np.array([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 10]
    ])
    B = np.array([
        [1, 0, 1],
        [0, 2, 0],
        [0, 0, 3]
    ])

    # Check if matrices are square and of the same dimension
    if A.shape != B.shape or A.shape[0] != A.shape[1]:
        print("Error: Matrices A and B must be square and of the same size.")
        return

    n = A.shape[0]

    # 1. Compute singular values of A and B.
    # np.linalg.svd returns singular values in descending order.
    a_vals = np.linalg.svd(A, compute_uv=False)
    b_vals = np.linalg.svd(B, compute_uv=False)

    # 2. Calculate the sign of det(A)*det(B)
    det_A = np.linalg.det(A)
    det_B = np.linalg.det(B)
    s = np.sign(det_A * det_B)

    # If a determinant is zero, its sign is 0.
    if s == 0:
        # Check if either matrix is singular, making a_n or b_n zero.
        # If a_n or b_n is approx zero, the last term vanishes anyway.
        # The sign function can be treated as 0 for this case.
        # np.sign(0) correctly returns 0.
        pass

    # 3. Apply the formula
    # min = - sum_{i=1}^{n-1} a_i*b_i - s * a_n*b_n
    sum_part = sum(a_vals[i] * b_vals[i] for i in range(n - 1))
    last_term = s * a_vals[n - 1] * b_vals[n - 1]
    min_value = -sum_part - last_term

    # 4. Print the results in a detailed way
    print("The minimum of the asymmetric Brockett cost function is given by the formula:")
    print(f"min = - (a_1*b_1 + ... + a_{n-1}*b_{n-1}) - s * a_n*b_n")
    print("where a_i and b_i are the singular values, and s = sign(det(A)*det(B)).")

    print("\nFor the given matrices A and B:")
    print("\nMatrix A:\n", A)
    print("\nMatrix B:\n", B)

    print("\nCalculated values:")
    # Using np.round for cleaner output
    print(f"  Singular values of A (a_i): {np.round(a_vals, 4)}")
    print(f"  Singular values of B (b_i): {np.round(b_vals, 4)}")
    print(f"  det(A) = {det_A:.4f}, det(B) = {det_B:.4f}")
    print(f"  s = sign({det_A:.4f} * {det_B:.4f}) = {s:.0f}")

    # Build and print the final equation string
    sum_expr_parts = []
    for i in range(n - 1):
        sum_expr_parts.append(f"{a_vals[i]:.4f} * {b_vals[i]:.4f}")
    
    print("\nFinal equation with numerical values:")
    sum_expression = " + ".join(sum_expr_parts)
    # The prompt requires outputting each number in the final equation.
    # The code below constructs the string showing the summation explicitly.
    equation_string = f"min = -({sum_expression}) - ({s:.0f} * {a_vals[n-1]:.4f} * {b_vals[n-1]:.4f})"
    print(equation_string)

    print("\nFinal Result:")
    print(f"min = {min_value}")

if __name__ == '__main__':
    solve_brockett_minimum()
