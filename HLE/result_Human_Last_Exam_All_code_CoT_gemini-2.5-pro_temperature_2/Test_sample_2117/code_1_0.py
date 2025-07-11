import numpy as np

def calculate_E(matrix):
    """Calculates the average eigenvalue gap of a matrix."""
    # Use eigvalsh for real symmetric matrices
    eigs = np.linalg.eigvalsh(matrix)
    eigs.sort()
    gaps = np.diff(eigs)
    # The definition is for the absolute difference
    return np.sum(np.abs(gaps)) / (len(eigs) - 1)

def calculate_S(matrix):
    """Calculates the mean square of the singular values of a matrix."""
    # SVD computes singular values, which are non-negative.
    s_vals_sq = np.linalg.svd(matrix, compute_uv=False)**2
    return np.mean(s_vals_sq)

def solve_and_verify(n):
    """
    Numerically verifies the derived formulas for a given integer n.
    """
    print(f"--- Verification for n = {n} ---")
    m = n + 2
    
    # 1. Cayley-Menger Matrix C
    C = np.ones((m, m)) - np.eye(m)
    
    # 2. H is unitarily similar to C, so E_H = E_C, S_H = S_C
    E_H_val = calculate_E(C)
    S_H_val = calculate_S(C)

    E_H_formula = (n + 2) / (n + 1)
    S_H_formula = float(n + 1)

    print(f"E_H (numerical from C) = {E_H_val:.6f}, E_H (formula: (n+2)/(n+1)) = {E_H_formula:.6f}")
    print(f"S_H (numerical from C) = {S_H_val:.6f}, S_H (formula: n+1)       = {S_H_formula:.6f}")
    
    # 3. Construct the Householder transformation matrix P
    # The transformation P that tridiagonalizes C is a single Householder matrix.
    P_hat = np.eye(m - 1)
    # Define vector for reflection
    x = np.ones(m - 1)
    # Build Householder reflection for x
    alpha = -np.sign(x[0]) * np.linalg.norm(x)
    u = x.copy()
    u[0] -= alpha
    v = u / np.linalg.norm(u)
    P_hat = np.eye(m - 1) - 2 * np.outer(v, v)

    # Embed P_hat into the full P matrix
    P = np.eye(m)
    P[1:, 1:] = P_hat

    # 4. Calculate E_P and S_P
    E_P_val = calculate_E(P)
    S_P_val = calculate_S(P)

    E_P_formula = 2 / (n + 1)
    S_P_formula = 1.0
    
    print(f"E_P (numerical)         = {E_P_val:.6f}, E_P (formula: 2/(n+1))   = {E_P_formula:.6f}")
    print(f"S_P (numerical)         = {S_P_val:.6f}, S_P (formula: 1)         = {S_P_formula:.6f}")

    # 5. Calculate the total product
    product_numerical = E_P_val * E_H_val * S_P_val * S_H_val
    product_formula = E_P_formula * E_H_formula * S_P_formula * S_H_formula
    print(f"\nProduct (numerical) = {product_numerical:.6f}")
    print(f"Product (formula)   = {product_formula:.6f}\n")


# --- Main Calculation ---
# First, verify the derived formulas for a sample value of n.
solve_and_verify(n=3)

# The product expression X(n) = 2(n+2)/(n+1) is a strictly decreasing function for n>=1.
# Therefore, its least upper bound (supremum) is its value at n=1.
n_for_lub = 1
numerator = 2 * (n_for_lub + 2)
denominator = n_for_lub + 1
lub = numerator / denominator

print("-" * 35)
print("The least upper bound of the product is found by maximizing the expression:")
print("f(n) = 2 * (n + 2) / (n + 1)")
print("This function decreases as n increases, so the maximum is at n = 1.")
print("\nFinal Calculation:")
print(f"For n = {n_for_lub}, the product is:")
print(f"  2 * ({n_for_lub} + 2) / ({n_for_lub} + 1) = (2 * 3) / 2 = {int(lub)}")
print(f"\nThe least upper bound of the product E_P * E_H * S_P * S_H is {int(lub)}.")
