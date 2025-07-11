import numpy as np

def solve_kk_eigenvalues():
    """
    This function solves the eigenvalue problem for the Kaluza-Klein modes
    of a spin-2 graviton in the specified 5D warped geometry.
    It counts the number of mass-squared eigenvalues below a given threshold.
    """
    
    # Parameters from the problem statement
    # A(x) = C1*sin(x) + C2*cos(x)
    C1 = 1.0
    C2 = 4.0
    # The threshold for the eigenvalues m^2
    THRESHOLD = 14.0

    # Numerical parameters
    # N is the number of grid points for discretization. A larger N gives higher accuracy.
    N = 500
    
    # Discretize the domain x in [0, 2*pi]
    # We use endpoint=False because of periodic boundary conditions x(0) = x(2*pi)
    x_grid = np.linspace(0, 2 * np.pi, N, endpoint=False)
    dx = x_grid[1] - x_grid[0]

    # Define the warp factor and its derivatives
    A = C1 * np.sin(x_grid) + C2 * np.cos(x_grid)
    A_prime = C1 * np.cos(x_grid) - C2 * np.sin(x_grid)
    A_double_prime = -C1 * np.sin(x_grid) - C2 * np.cos(x_grid)

    # Define the coefficients of the differential operator H:
    # H = p(x) * d^2/dx^2 + q(x) * d/dx + r(x)
    p = -np.exp(-2 * A)
    q = A_prime * np.exp(-2 * A)
    r = (1.5 * A_double_prime + 0.75 * A_prime**2) * np.exp(-2 * A)

    # Construct the matrix for the discretized operator H
    M = np.zeros((N, N), dtype=float)
    
    # Fill the matrix using central differences and periodic boundary conditions
    for i in range(N):
        # Indices for neighbors with periodic boundary conditions
        prev_i = (i - 1 + N) % N
        next_i = (i + 1) % N
        
        # Diagonal term
        M[i, i] = -2 * p[i] / (dx**2) + r[i]
        
        # Off-diagonal terms
        M[i, next_i] = p[i] / (dx**2) + q[i] / (2 * dx)
        M[i, prev_i] = p[i] / (dx**2) - q[i] / (2 * dx)

    # Compute the eigenvalues of the matrix M
    eigenvalues = np.linalg.eigvals(M)
    
    # Eigenvalues should be real, but numerical methods might introduce small imaginary parts
    real_eigenvalues = np.real(eigenvalues)

    # Filter for eigenvalues below the threshold and sort them
    eigenvalues_below_threshold = np.sort(real_eigenvalues[real_eigenvalues < THRESHOLD])
    
    # Count the number of such eigenvalues
    count = len(eigenvalues_below_threshold)

    # Output the results as requested. 
    # The "final equation" shows which numbers contribute to the final count.
    print(f"The mass-squared eigenvalues below {THRESHOLD} are:")
    # Create the equation string showing each value found.
    equation_parts = [f"{eig:.4f}" for eig in eigenvalues_below_threshold]
    if not equation_parts:
        print("None")
        final_equation = "0"
    else:
        # We output each eigenvalue that is part of the final count.
        # Let's show this as an equation summing up 1 for each found eigenvalue.
        sum_str = " + ".join(["1"] * count)
        final_equation = f"{sum_str} = {count}"
        print(np.round(eigenvalues_below_threshold, 4))

    print(f"\nThe final equation for the count is: {final_equation}")
    print(f"\nTotal count of eigenvalues below {THRESHOLD}: {count}")

    return count

if __name__ == '__main__':
    final_count = solve_kk_eigenvalues()
    # The final answer in the requested format
    # <<<final_count>>>
    # Since this is a self-contained script, we will print the final answer
    # to be easily extracted.

<<<solve_kk_eigenvalues()>>>