import numpy as np

def solve_phi_n(n=5):
    """
    Calculates the value of phi(n) based on the derived analytical solution
    and verifies it through numerical matrix calculations.

    The problem simplifies to calculating phi(n) = exp(tr(Proj(X_inv))),
    which further simplifies to exp(2n - 2/n). This function demonstrates
    the steps to reach this result numerically.
    """
    if n < 5:
        print("Error: n must be greater than or equal to 5.")
        return

    # Step 1: Define the vector u
    # The integrals in the problem evaluate to I1 = pi/2 and I2 = -pi/2.
    # [u]_i = -pi/2 for i odd, pi/2 for i even (using 1-based indexing).
    pi_half = np.pi / 2
    # Create u vector based on 1-based indexing for i from 1 to n
    u = np.array([pi_half * (-1)**i for i in range(1, n + 1)])

    # Step 2: Define the matrix Y = X^(-1)
    # The matrix X is defined such that its inverse, Y, is a simple
    # symmetric tridiagonal matrix with 2 on the diagonal and 1 on
    # the super- and sub-diagonals.
    Y = np.zeros((n, n))
    Y.flat[::n+1] = 2  # Set diagonal to 2
    Y.flat[1::n+1] = 1  # Set superdiagonal to 1
    Y.flat[n::n+1] = 1  # Set subdiagonal to 1

    # Step 3: Calculate the components for the trace of the projected matrix Z = Proj(Y)
    # The trace is given by tr(Z) = tr(Y) - tr(W), where W is the projection
    # of Y onto the orthogonal complement of the tangent space.
    # tr(W) can be calculated as (u^T * Y * u) / ||u||^2.

    # Calculate tr(Y)
    tr_Y = np.trace(Y)

    # Calculate the quadratic form u^T * Y * u
    uT_Y_u = u @ Y @ u

    # Calculate the squared norm of u, ||u||^2
    norm_u_sq = u @ u

    # Calculate tr(W)
    tr_W = uT_Y_u / norm_u_sq

    # Calculate tr(Z)
    tr_Z = tr_Y - tr_W

    # Step 4: Calculate the final determinant value
    # phi(n) = det(Expm(Z)) = exp(tr(Z))
    phi_n = np.exp(tr_Z)

    # Output the results, showing each number in the final equation.
    print(f"Calculation for n = {n}")
    print("-" * 30)
    print("The final value is derived from the formula: phi(n) = exp(tr(Z))")
    print("where tr(Z) = tr(Y) - tr(W), and Y = X_inv.")
    print("\nThe components of the trace calculation are:")
    print(f"tr(Y) = 2*n = {tr_Y:.4f}")
    # The expression for tr(W) involves several components:
    print(f"tr(W) is calculated from (u^T * Y * u) / ||u||^2:")
    print(f"  - The quadratic form u^T * Y * u = {uT_Y_u:.4f}")
    print(f"  - The squared norm ||u||^2 = {norm_u_sq:.4f}")
    print(f"  - Resulting tr(W) = 2/n = {tr_W:.4f}")
    print(f"\nTrace of the projected matrix Z:")
    print(f"tr(Z) = tr(Y) - tr(W) = {tr_Y:.1f} - {tr_W:.1f} = {tr_Z:.4f}")
    print("\nFinal determinant value:")
    print(f"phi({n}) = exp({tr_Z:.4f}) = {phi_n}")

if __name__ == '__main__':
    # The problem is stated for n >= 5. We use n=5 as an example.
    solve_phi_n(n=5)