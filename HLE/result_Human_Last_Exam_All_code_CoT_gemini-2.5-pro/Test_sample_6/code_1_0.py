import numpy as np
from scipy.linalg import eigh

def solve_kk_masses():
    """
    This function calculates the Kaluza-Klein mass eigenvalues for a given
    5D gravitational theory and counts how many are below a threshold.
    """
    # 1. Setup the problem parameters
    N = 500  # Number of grid points for discretization
    L = 2 * np.pi  # Length of the extra dimension
    h = L / N  # Grid spacing
    threshold = 14

    # Create the grid for the extra dimension coordinate x
    # Main grid for psi and w(x)
    x = np.linspace(0, L, N, endpoint=False)
    # Staggered grid for p(x) to improve accuracy of the finite difference method
    x_staggered = x + h / 2

    # 2. Define the functions from the Sturm-Liouville problem
    # The warp factor A(x)
    def A(t):
        return np.sin(t) + 4 * np.cos(t)

    # The functions in the eigenvalue problem: (p(x)psi')' + m^2*w(x)*psi = 0
    # For this problem, p(x) = w(x) = exp(3*A(x))
    def p_and_w(t):
        return np.exp(3 * A(t))

    # Calculate p on the staggered grid and w on the main grid
    p_vals_staggered = p_and_w(x_staggered)
    w_vals = p_and_w(x)

    # 3. Construct the matrices K and W for the generalized eigenvalue problem K*v = lambda*W*v
    K = np.zeros((N, N))
    h2 = h**2

    # Fill the matrix K representing the kinetic operator -(p(x)psi')'
    # This uses central differences on a staggered grid for stability and accuracy.
    for i in range(N):
        # Diagonal elements
        p_plus_half = p_vals_staggered[i]
        p_minus_half = p_vals_staggered[i - 1]  # Periodic boundary handled by numpy's negative index
        K[i, i] = (p_plus_half + p_minus_half) / h2

        # Off-diagonal elements (exploiting K's symmetry)
        K[i, (i + 1) % N] -= p_plus_half / h2
        K[(i + 1) % N, i] -= p_plus_half / h2
    
    # The matrix W is diagonal with the values of w(x)
    W = np.diag(w_vals)

    # 4. Solve the generalized eigenvalue problem
    # The function eigh is for symmetric/Hermitian matrices and is numerically stable.
    # The returned eigenvalues correspond to m^2.
    eigenvalues, eigenvectors = eigh(K, W)

    # 5. Count the eigenvalues below the threshold and print the results
    eigenvalues_below_threshold = []
    count = 0
    for ev in sorted(eigenvalues):
        if ev < threshold:
            count += 1
            eigenvalues_below_threshold.append(ev)

    print("The eigenvalues m^2 below 14 are:")
    for ev in eigenvalues_below_threshold:
        print(f"{ev:.4f}")
    
    print(f"\nThe number of eigenvalues below {threshold} is: {count}")

    # Return the final count for the answer format
    return count

if __name__ == '__main__':
    final_count = solve_kk_masses()
    # The final answer format is not printed here as it's part of the assistant's response wrapper.
    # To run this script standalone, you can uncomment the next line.
    # print(f"<<<{final_count}>>>")

# Execute the function to produce the output
final_count = solve_kk_masses()
print(f"\n<<<>>>\n{final_count}")