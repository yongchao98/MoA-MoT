import numpy as np

def compute_ngd_update_efficiently(d, n, alpha, X, g_vec):
    """
    Computes the NGD update direction p = (F + alpha*I)^-1 * g efficiently.

    Args:
        d (int): Dimension of the weight matrix.
        n (int): Number of samples.
        alpha (float): Damping factor.
        X (np.ndarray): Input data matrix of shape (d, n).
        g_vec (np.ndarray): Gradient vector of shape (d*d,).

    Returns:
        np.ndarray: The update direction vector p of shape (d*d,).
    """
    # The cost of this algorithm is O(d^2 * n), which is much better than
    # the naive O(d^6).

    # Reshape the gradient vector g into a d x d matrix G for easier handling.
    G = g_vec.reshape(d, d)

    # p = (1/alpha)g - (1/alpha^2) * (X kron I) * M^-1 * (X.T kron I) * g
    # where M = I_nd + (1/alpha) * (X.T*X kron I_d)

    # 1. Compute C_n = X.T @ X. Cost: O(n^2 * d)
    C_n = X.T @ X

    # 2. Eigendecomposition of C_n = Q*S*Q.T. Cost: O(n^3)
    s_diag, Q = np.linalg.eigh(C_n)

    # 3. Compute the term (X.T kron I_d) * g. In matrix form, this is G @ X.
    # Cost: O(d^2 * n)
    U_mat = G @ X

    # 4. We need to compute z = M^-1 * u, where u = vec(U_mat).
    # z = (Q kron I) * D_S * (Q.T kron I) * u
    # where D_S is a diagonal matrix.
    # This can be done with matrix operations.
    # z_intermediate = U_mat @ Q. Cost: O(d * n^2)
    z_intermediate = U_mat @ Q

    # Scale columns by inverse eigenvalues of M. Cost: O(d * n)
    scaling_factors = 1.0 / (1.0 + s_diag / alpha)
    z_scaled = z_intermediate * scaling_factors # Broadcasting

    # z_mat = z_scaled @ Q.T. Cost: O(d * n^2)
    Z_mat = z_scaled @ Q.T

    # 5. Compute w = (X kron I_d) * z. In matrix form, this is Z_mat @ X.T
    # Cost: O(d^2 * n)
    W_mat = Z_mat @ X.T

    # 6. Combine terms to get the final update matrix P. Cost: O(d^2)
    P_mat = (1.0 / alpha) * G - (1.0 / alpha**2) * W_mat

    return P_mat.flatten()

if __name__ == '__main__':
    # --- Problem Setup ---
    # Set dimensions d and number of samples n, with n < d
    d = 100
    n = 20
    alpha = 0.1 # Damping factor

    # --- Generate Random Data for Demonstration ---
    # Generate random input data matrix X
    X = np.random.randn(d, n)
    # Generate a random gradient vector g
    g_vec = np.random.randn(d * d)

    # --- Compute the Update ---
    # This call executes the efficient O(d^2 * n) algorithm
    p_efficient = compute_ngd_update_efficiently(d, n, alpha, X, g_vec)

    # --- Output the Explanation and Result ---
    print("--- Natural Gradient Descent Curvature Cost Analysis ---")
    print(f"For a network layer of size d x d and n samples (d={d}, n={n}):\n")
    print("The Fisher Information Matrix F is d^2 x d^2, but has a low-rank structure.")
    print("A direct inversion of (F + alpha*I) would cost O(d^6).\n")
    print("By exploiting the structure F = (X*X^T) kron I_d using the Woodbury identity,")
    print("the cost can be significantly reduced.\n")

    print("--- Minimum Curvature Cost Equation ---")
    final_equation = "O(d^exponent_d * n^exponent_n)"
    exponent_d = 2
    exponent_n = 1
    print(f"The minimum achievable curvature cost is described by the equation: {final_equation}")
    print(f"The numbers in this final equation are:")
    print(f"  d = {d}")
    print(f"  exponent_d = {exponent_d}")
    print(f"  n = {n}")
    print(f"  exponent_n = {exponent_n}\n")

    print(f"The complexity is O(d^{exponent_d} * n^{exponent_n}).\n")

    print("--- Result of the Computation ---")
    print("The NGD update direction vector p has been computed.")
    print(f"Shape of p: {p_efficient.shape}")
    print(f"Norm of p: {np.linalg.norm(p_efficient):.4f}")
