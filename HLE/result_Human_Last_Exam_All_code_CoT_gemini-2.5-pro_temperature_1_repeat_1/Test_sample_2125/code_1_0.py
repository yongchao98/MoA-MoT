import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import root_scalar

def solve_for_E2(alpha):
    """
    Solves the time-independent Schrödinger equation for the given potential
    and returns the second excited state energy (E2).

    Args:
        alpha (float): The parameter alpha in the potential.

    Returns:
        float: The energy of the second excited state.
    """
    # Define the spatial grid
    N = 2000  # Number of grid points
    L = 10.0  # Grid extends from -L to L
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]

    # Define the potential V(x)
    V = 0.5 * x**6 - alpha * x**4 + 0.5 * (alpha**2 - 7) * x**2

    # Construct the Hamiltonian matrix using the finite difference method.
    # The Hamiltonian is a tridiagonal matrix.
    diagonal = 1/dx**2 + V
    off_diagonal = -1/(2*dx**2) * np.ones(N - 1)

    # Find the eigenvalues. We are interested in the low-lying states.
    # eigh_tridiagonal is efficient for this.
    # We select the first 3 eigenvalues (indices 0, 1, 2).
    eigenvalues, _ = eigh_tridiagonal(diagonal, off_diagonal, select='i', select_range=(0, 2))

    # The second excited state energy E2 corresponds to the third eigenvalue (index 2).
    return eigenvalues[2]

# We need to find the value of alpha for which E2(alpha) = 0.
# We can observe that for alpha > sqrt(7), the potential wells deepen, and
# the energy levels decrease. E2 will eventually cross zero.
# Let's find a bracket for the root.
# For alpha = 2.7, E2 is positive.
# For alpha = 3.0, E2 is negative.
# So the root lies between 2.7 and 3.0.
try:
    solution = root_scalar(solve_for_E2, bracket=[2.7, 3.0], method='brentq')
    alpha_0 = solution.root
    e2_at_alpha_0 = solve_for_E2(alpha_0)

    # Print the result in the requested format.
    # The final equation is E_2(alpha_0) = 0, which makes F(alpha_0) = 0.
    print(f"The largest value alpha_0 is found where the second excited state energy E_2 is zero.")
    print(f"The condition is E_2({alpha_0:.4f}) = {e2_at_alpha_0:.4f}")
    print(f"\nThus, the largest value is alpha_0 = {alpha_0:.4f}")
    print(f'<<<{alpha_0:.4f}>>>')

except (ValueError, RuntimeError) as e:
    print(f"An error occurred during root finding: {e}")
    print("Could not determine the value of alpha_0.")
