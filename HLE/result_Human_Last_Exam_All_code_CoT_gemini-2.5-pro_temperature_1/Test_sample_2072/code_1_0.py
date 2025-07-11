import math

def solve_phi(n):
    """
    Calculates the value of phi(n) based on the derived formula.

    The derivation proceeds as follows:
    1.  The problem statement implies that for M_ij = 1/n, M must be on the manifold, which necessitates I1 = I2.
        This makes the vectors u and v constant vectors (multiples of 1_n).
    2.  The tangent space is then {Z in SYM(n) | Z * 1_n = 0}.
    3.  The value to compute is phi(n) = det(Expm(Proj_M(X_inv))) = exp(Tr(Proj_M(X_inv))).
    4.  Tr(Proj_M(Y)) for Y in SYM(n) is Tr(Y) - (1/n) * 1_n^T * Y * 1_n.
    5.  The matrix X is found to be X = diag(h) * T_inv * diag(h), where T is the tridiagonal matrix [2, -1; -1, 2, ...].
    6.  This gives X_inv = diag(h) * T * diag(h), which is a tridiagonal matrix with 2 on the diagonal and 1 on the off-diagonals.
    7.  Tr(X_inv) = 2n.
    8.  1_n^T * X_inv * 1_n = 4n - 2.
    9.  Tr(Proj_M(X_inv)) = 2n - (4n - 2)/n = 2*(n-1)**2 / n.
    10. phi(n) = exp(2*(n-1)**2 / n).
    """
    if n < 5:
        raise ValueError("The problem is defined for n >= 5.")

    # Calculate the trace of the projection
    trace_proj = 2 * (n - 1)**2 / n

    # Calculate phi(n)
    phi_n = math.exp(trace_proj)

    # Print the final equation as requested
    print(f"For n = {n}:")
    print(f"Tr(Proj_M(X^-1)) = 2*({n}-1)^2 / {n} = {trace_proj}")
    print(f"phi({n}) = exp({trace_proj}) = {phi_n}")
    
    return phi_n

# The problem is stated for n >= 5. Let's calculate for n=5.
n_val = 5
final_value = solve_phi(n_val)
# The final answer is the numerical value for n=5.
# print(f"<<<{final_value}>>>") # This would be the submission format