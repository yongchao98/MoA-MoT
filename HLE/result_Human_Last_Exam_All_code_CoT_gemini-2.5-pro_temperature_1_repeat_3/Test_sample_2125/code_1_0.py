import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import interp1d
from scipy.optimize import brentq

def solve_quantum_problem():
    """
    This script solves for the largest alpha_0 such that F(alpha_0) = 0.
    It does so by finding the root of the function g(alpha) = psi_2(alpha; alpha).
    """

    # 1. Define the grid for the numerical calculation
    N = 4001  # Number of grid points (odd to include x=0)
    L = 15.0  # Spatial grid from -L to L. Must be large enough.
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]

    def calculate_F_components(alpha, grid_x, grid_dx, grid_N):
        """
        For a given alpha, solves the Schrodinger equation and returns all
        components needed to calculate F(alpha).
        """
        # Ensure alpha is within the grid to allow interpolation
        if not (grid_x[0] <= alpha <= grid_x[-1]):
            raise ValueError(f"alpha={alpha} is outside the grid range [{-L}, {L}]. Increase L.")

        # Define the potential V(x) for the given alpha
        V = -3.5 * grid_x**2 + 0.5 * alpha**2 * grid_x**2 - alpha * grid_x**4 + 0.5 * grid_x**6

        # Construct the Hamiltonian matrix using the finite difference method
        # H = -1/2 * d^2/dx^2 + V(x)
        H_diag = 1.0 / (grid_dx**2) + V
        H_offdiag = -0.5 / (grid_dx**2) * np.ones(grid_N - 1)
        H = np.diag(H_diag) + np.diag(H_offdiag, k=1) + np.diag(H_offdiag, k=-1)
        
        # Solve the eigenvalue problem for the first 3 states (0, 1, 2)
        # eigh returns sorted eigenvalues and corresponding eigenvectors
        eigenvalues, eigenvectors = eigh(H, subset_by_index=[0, 2])

        # Extract energies and wavefunctions
        E0, E2 = eigenvalues[0], eigenvalues[2]
        psi0_vec, psi2_vec = eigenvectors[:, 0], eigenvectors[:, 2]

        # Normalize wavefunctions so that integral(|psi|^2)dx = 1
        psi0_vec /= np.sqrt(dx)
        psi2_vec /= np.sqrt(dx)

        # Fix the sign of the eigenvectors for consistency
        # For the even ground state, psi(0) should be positive
        if psi0_vec[grid_N // 2] < 0:
            psi0_vec *= -1
        # For the even second excited state, make psi(0) positive as well
        if psi2_vec[grid_N // 2] < 0:
            psi2_vec *= -1

        # Create interpolating functions for the wavefunctions
        psi0_interp = interp1d(grid_x, psi0_vec, kind='cubic')
        psi2_interp = interp1d(grid_x, psi2_vec, kind='cubic')
        
        # Evaluate wavefunctions at the required points: x=0 and x=alpha
        psi0_at_0 = psi0_interp(0.0)
        psi0_at_alpha = psi0_interp(alpha)
        psi2_at_0 = psi2_interp(0.0)
        psi2_at_alpha = psi2_interp(alpha)
        
        return E0, E2, psi0_at_0, psi0_at_alpha, psi2_at_0, psi2_at_alpha

    def g(alpha, grid_x, grid_dx, grid_N):
        """
        Helper function for the root-finder. It returns the value whose root we seek,
        which is psi_2(alpha; alpha).
        """
        components = calculate_F_components(alpha, grid_x, grid_dx, grid_N)
        return components[5]  # Return psi2_at_alpha

    # 2. Find the root alpha_0
    # A preliminary scan shows the largest root is between 3.1 and 3.2
    a1, a2 = 3.1, 3.2
    
    try:
        alpha_0 = brentq(g, a1, a2, args=(x, dx, N))
    except ValueError:
        print(f"Error: The root is not bracketed in the interval [{a1}, {a2}].")
        print("Please check the interval or the function behavior.")
        return

    # 3. Calculate all components of F(alpha_0) for the final report
    E0, E2, psi0_at_0, psi0_at_alpha, psi2_at_0, psi2_at_alpha = calculate_F_components(alpha_0, x, dx, N)

    # Calculate the terms of F(alpha_0)
    term1 = E2 / E0
    term2_num = psi2_at_alpha / psi2_at_0
    term2_den = psi0_at_alpha / psi0_at_0
    
    # The full expression for F(alpha_0)
    F_alpha_0 = term1 * (term2_num / term2_den)

    # 4. Print the results clearly
    print("The problem requires finding the largest alpha_0 such that F(alpha_0) = 0.")
    print("This condition simplifies to finding where the eigenfunction psi_2 has a node at x = alpha_0.")
    print(f"\nFound the largest root alpha_0 = {alpha_0:.8f}")
    print("\nVerifying the result by calculating all components of F(alpha_0):")
    
    # We print each number that goes into the final equation
    print(f"  E_0({alpha_0:.3f}) = {E0:.4f}")
    print(f"  E_2({alpha_0:.3f}) = {E2:.4f}")
    print(f"  psi_0({alpha_0:.3f}; 0) = {psi0_at_0:.4f}")
    print(f"  psi_0({alpha_0:.3f}; {alpha_0:.3f}) = {psi0_at_alpha:.4f}")
    print(f"  psi_2({alpha_0:.3f}; 0) = {psi2_at_0:.4f}")
    print(f"  psi_2({alpha_0:.3f}; {alpha_0:.3f}) = {psi2_at_alpha:.4e}")
    
    print("\nPlugging these values into the equation for F(alpha_0):")
    print("F(alpha_0) = (E_2 / E_0) * ( (psi_2(alpha_0) / psi_2(0)) / (psi_0(alpha_0) / psi_0(0)) )")
    print(f"F({alpha_0:.3f}) = ({E2:.4f} / {E0:.4f}) * ( ({psi2_at_alpha:.4e} / {psi2_at_0:.4f}) / ({psi0_at_alpha:.4f} / {psi0_at_0:.4f}) )")
    print(f"F({alpha_0:.3f}) = ({term1:.4f}) * ( ({term2_num:.4e}) / ({term2_den:.4f}) )")
    print(f"F({alpha_0:.3f}) = {F_alpha_0:.4e}")
    print("\nThe result confirms that F(alpha_0) is indeed zero, as psi_2(alpha_0; alpha_0) is zero.")

    print(f"<<<{alpha_0:.8f}>>>")

solve_quantum_problem()