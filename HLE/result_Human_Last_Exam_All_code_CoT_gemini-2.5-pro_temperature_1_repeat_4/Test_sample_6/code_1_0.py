import numpy as np
from scipy.linalg import eig

def solve_and_count_eigenvalues():
    """
    Solves the Kaluza-Klein mass eigenvalue problem numerically and counts
    the number of eigenvalues below a specified threshold.
    """
    # 1. Define problem parameters
    N = 1000  # Number of grid points for discretization
    threshold = 14.0

    # 2. Set up the grid and define the relevant functions from the theory
    x, h = np.linspace(0, 2 * np.pi, N, endpoint=False, retstep=True)

    def A(x_val):
        """The warp factor function."""
        return np.sin(x_val) + 4 * np.cos(x_val)

    def p_func(x_val):
        """The function p(x) in the Sturm-Liouville operator."""
        return np.exp(3 * A(x_val))

    def r_func(x_val):
        """The weight function r(x) in the Sturm-Liouville problem."""
        return np.exp(A(x_val))

    # 3. Discretize the operator and construct the matrices M and R
    # We evaluate p(x) at points halfway between the grid points for better accuracy.
    p_half = p_func(x + h / 2)
    r_vals = r_func(x)

    # Construct the "stiffness" matrix M using finite differences
    # for the operator L = -(p(x)u')'
    diag_main = (p_half + np.roll(p_half, 1)) / h**2
    diag_off = -p_half / h**2

    M = np.diag(diag_main) + np.diag(diag_off[:-1], 1) + np.diag(diag_off[:-1], -1)
    
    # Apply periodic boundary conditions
    M[0, N - 1] = diag_off[-1]
    M[N - 1, 0] = diag_off[-1]

    # Construct the "mass" or "weight" matrix R
    R = np.diag(r_vals)

    # 4. Solve the generalized eigenvalue problem M*v = lambda*R*v
    eigenvalues, _ = eig(M, b=R)
    
    # The eigenvalues m^2 should be real. We take the real part to discard
    # small imaginary parts arising from numerical inaccuracies.
    m_squared = np.real(eigenvalues)

    # 5. Count the eigenvalues below the threshold
    count = np.sum(m_squared < threshold)

    print(f"The differential equation for the Kaluza-Klein modes psi(x) is of the form:")
    print("-(p(x) * psi'(x))' = m^2 * r(x) * psi(x)")
    print(f"where p(x)=exp(3*A(x)), r(x)=exp(A(x)), and A(x)=sin(x)+4*cos(x).")
    print(f"We are counting the number of eigenvalues m^2 below the threshold of {threshold}.")
    print("\n--- Result ---")
    print(f"The number of eigenvalues below {threshold} is: {count}")
    
    # As requested by the prompt format, printing a final equation is not directly
    # applicable, so we print the final count as the result.
    # For clarity, we also list the found eigenvalues below the threshold.
    m_squared.sort()
    found_eigenvalues = m_squared[m_squared < threshold]
    print(f"\nThe {len(found_eigenvalues)} eigenvalues m^2 found below {threshold} are:")
    # Rounding for cleaner output
    print(np.round(found_eigenvalues, 4))


solve_and_count_eigenvalues()
<<<33>>>