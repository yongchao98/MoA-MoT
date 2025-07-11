import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

def solve_for_largest_alpha():
    """
    Finds the largest value of alpha, alpha_0, for which F(alpha) = 0.
    This is achieved by finding the roots of psi_2(alpha, alpha) = 0.
    """

    def get_psi2_at_alpha(alpha):
        """
        Calculates the value of the 2nd excited state eigenfunction psi_2 at x=alpha.
        
        Args:
            alpha (float): The parameter of the potential.
            
        Returns:
            float: The value of psi_2(alpha, alpha).
        """
        # Set up a sufficiently large and dense grid for numerical accuracy.
        # The grid size L is made adaptive to ensure alpha is always within the boundaries.
        N = 4001
        L = max(15.0, alpha * 2.0)
        x_grid = np.linspace(-L, L, N)
        dx = x_grid[1] - x_grid[0]

        # Define the sextic anharmonic potential V(x)
        potential = (0.5 * alpha**2 - 3.5) * x_grid**2 - alpha * x_grid**4 + 0.5 * x_grid**6

        # Construct the Hamiltonian matrix using finite differences
        H_diag = 1/dx**2 + potential
        H_offdiag = -0.5/dx**2 * np.ones(N - 1)
        hamiltonian = np.diag(H_diag) + np.diag(H_offdiag, k=1) + np.diag(H_offdiag, k=-1)

        # Solve the eigenvalue problem for the first 3 states (n=0, 1, 2)
        # eigh is used for symmetric matrices; it's fast and sorts eigenvalues.
        eigenvalues, eigenvectors = eigh(hamiltonian, subset_by_index=[0, 2])
        
        # The 2nd excited state (psi_2) corresponds to the 3rd eigenvector (index 2)
        psi2_vec = eigenvectors[:, 2]

        # For consistency across different alpha values, we enforce a sign convention.
        # Since psi_2 is an even function, we can demand its value at x=0 to be positive.
        psi2_at_zero = np.interp(0, x_grid, psi2_vec)
        if psi2_at_zero < 0:
            psi2_vec = -psi2_vec

        # Use interpolation to find the value of psi_2 accurately at x = alpha
        psi2_at_alpha = np.interp(alpha, x_grid, psi2_vec)

        return psi2_at_alpha

    #--- Root Finding Procedure ---

    # 1. Perform a coarse scan to find intervals where psi_2(alpha, alpha) changes sign
    scan_min = 0.1
    scan_max = 6.0
    num_points = 200
    alphas_to_scan = np.linspace(scan_min, scan_max, num_points)
    
    # We evaluate our target function on the scanned alpha values
    f_values = np.array([get_psi2_at_alpha(a) for a in alphas_to_scan])

    # 2. Identify the intervals where a root may exist (sign change)
    root_intervals = []
    for i in range(num_points - 1):
        if np.sign(f_values[i]) != np.sign(f_values[i+1]):
            root_intervals.append((alphas_to_scan[i], alphas_to_scan[i+1]))

    # 3. Use a precise numerical solver (brentq) to find the roots in these intervals
    roots = []
    for a, b in root_intervals:
        try:
            root = brentq(get_psi2_at_alpha, a, b)
            roots.append(root)
        except (ValueError, RuntimeError):
            # This may happen if the signs are the same due to numerical noise, or if the function is not well-behaved
            # In our case, it's safe to just skip that interval.
            pass

    # 4. Find the largest root, as requested by the problem
    if roots:
        largest_root = max(roots)
        # We present the final answer clearly. The final number in the equation is the value of alpha.
        print(f"The largest value alpha_0 such that F(alpha_0) = 0 is found to be:")
        print(f"F({largest_root}) = 0")
    else:
        print("No roots were found in the scanned range.")

if __name__ == '__main__':
    solve_for_largest_alpha()
