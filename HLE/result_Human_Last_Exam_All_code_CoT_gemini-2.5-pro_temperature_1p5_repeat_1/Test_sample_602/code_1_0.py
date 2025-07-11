import numpy as np

def solve_for_l_n():
    """
    This script calculates the value of l(n) for a given n, following the steps
    outlined in the problem description. It then prints a detailed breakdown
    using a derived analytical formula.
    """
    # Per the problem, n must be an integer >= 5. We use n=5 as an example.
    n = 5

    print(f"### Starting Calculation for n = {n} ###")

    # --- Part 1: Numerical Calculation via Matrix Operations ---
    # This part follows the definition of l(n) step-by-step.

    # Step 1: Construct the matrix M_n
    # M_ij = alpha if i=j, beta otherwise
    alpha = np.sqrt(1 - (n - 1) / n**2)
    beta = 1 / n
    M = np.full((n, n), beta, dtype=np.float64)
    np.fill_diagonal(M, alpha)

    # Step 2: Construct the matrix P_n
    # P_ij = (-1)^(i+j) * (min(i,j) - ij/(n+1))
    P = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            # Convert 0-based index to 1-based for formula
            i1, j1 = i + 1, j + 1
            term = min(i1, j1) - (i1 * j1) / (n + 1)
            P[i, j] = ((-1)**(i1 + j1)) * term

    # Step 3: Compute Q = f^(3)(P_n), which is the inverse of P_n
    # The matrix Q is f^(3)(P_n). Since P_n is invertible, Q = P_n^-1.
    try:
        Q = np.linalg.inv(P)
    except np.linalg.LinAlgError:
        print("P_n is singular, falling back to pseudoinverse.")
        Q = np.linalg.pinv(P)

    # Step 4: Compute R = f_M^(2)(Q), the projection of Q onto the tangent space at M.
    # The formula is R = Q - M * mdiag(M.T @ Q)
    # The diagonal of (M.T @ Q) can be computed efficiently.
    D_diag = np.einsum('ji,jk->k', M, Q)
    D = np.diag(D_diag)
    R = Q - M @ D

    # Step 5: Compute l(n) = f^(1)(R), the sum of the first and last rows of R.
    l_n_numerical = np.sum(R[0, :]) + np.sum(R[n - 1, :])


    # --- Part 2: Printout Using the Analytical Formula ---
    # To satisfy the "output each number in the final equation" requirement,
    # we use the derived analytical formula for l(n).
    
    print("\n### Detailed Breakdown of the Final Equation ###")
    print("The exact value of l(n) can be found using the analytical formula:")
    print("l(n) = (2 * (n^2 + 1 - (2*n - 1) * sqrt(n^2 - n + 1))) / n^2\n")

    print(f"For n = {n}, let's compute the components of the equation:")

    comp_n_sq = n**2
    print(f"n^2 = {n}^2 = {comp_n_sq}")

    comp1 = n**2 + 1
    print(f"n^2 + 1 = {comp_n_sq} + 1 = {comp1}")

    comp2 = 2*n - 1
    print(f"2*n - 1 = 2*{n} - 1 = {comp2}")

    comp3 = n**2 - n + 1
    print(f"n^2 - n + 1 = {comp_n_sq} - {n} + 1 = {comp3}")
    
    # Assemble the final equation string with the computed numbers.
    num_part1_val = 2 * comp1
    num_part2_val = -2 * comp2
    den_val = n**2

    print(f"\nSubstituting these values results in the final equation:")
    print(f"l({n}) = (2 * ({comp1} - {comp2} * sqrt({comp3}))) / {den_val}")
    print(f"     = ({num_part1_val} + ({num_part2_val}) * sqrt({comp3})) / {den_val}")
    
    # Calculate the final value from the analytical expression
    l_n_analytical = (num_part1_val + num_part2_val * np.sqrt(comp3)) / den_val
    
    print("\n--- Final Answer ---")
    print(f"The numerical result from matrix operations is: l({n}) = {l_n_numerical}")
    print(f"The exact value derived from the formula is: l({n}) = {l_n_analytical}")
    print(f"(The values should be nearly identical, differing only by floating-point precision.)")

if __name__ == '__main__':
    solve_for_l_n()