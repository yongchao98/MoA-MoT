import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

def get_E2(alpha: float) -> float:
    """
    Solves the 1D time-independent Schrodinger equation for a given alpha
    and returns the energy of the second excited state (E_2).

    Args:
        alpha: The positive parameter in the potential.

    Returns:
        The energy E_2.
    """
    # Set up the spatial grid. A larger grid is needed for larger alpha,
    # but L=15 is sufficient for the region of interest.
    # A high number of points ensures accuracy.
    N = 2001
    L = 15.0
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]

    # Define the potential V(x) for the given alpha
    # V(x) = (-7/2)x^2 + (1/2)alpha^2 x^2 - alpha x^4 + (1/2)x^6
    V = (-3.5 * x**2) + (0.5 * alpha**2 * x**2) - (alpha * x**4) + (0.5 * x**6)

    # Construct the Hamiltonian matrix using the finite difference method.
    # The second derivative is approximated as (psi_{i+1} - 2*psi_i + psi_{i-1}) / dx^2.
    # The kinetic energy part is -1/2 * d^2/dx^2.
    main_diag = 1.0 / dx**2 + V
    off_diag = -0.5 / dx**2 * np.ones(N - 1)
    
    # Use eigh to find eigenvalues and eigenvectors of the symmetric Hamiltonian matrix.
    # eigh is efficient and returns sorted eigenvalues.
    energies = eigh_tridiagonal(main_diag, off_diag, eigvals_only=True)
    
    # The energy levels are sorted in ascending order.
    # E_0 is energies[0], E_1 is energies[1], and E_2 is energies[2].
    E2 = energies[2]
    
    return E2

def solve():
    """
    Finds the value of alpha_0 for which F(alpha_0) = 0 by finding the root
    of E_2(alpha) = 0.
    """
    # The condition F(alpha) = 0 is satisfied if E_2(alpha) = 0.
    # We will find the value of alpha for which this is true.
    
    # Based on a preliminary scan, the root for E_2(alpha) = 0 lies
    # between alpha = 5.0 and alpha = 6.0.
    # E_2(5.0) is positive and E_2(6.0) is negative.
    lower_bound = 5.0
    upper_bound = 6.0

    # Use the Brent method, a reliable root-finding algorithm.
    try:
        alpha_0 = brentq(get_E2, lower_bound, upper_bound, xtol=1e-12, rtol=1e-12)
        print(f"The largest value alpha_0 such that F(alpha_0) = 0 is found by solving E_2(alpha_0) = 0.")
        # The final output needs to include each number of the equation.
        # Here we present the discovered value.
        print(f"alpha_0 = {alpha_0}")
    except ValueError:
        print("Failed to find a root in the given interval.")

# SciPy's eigh_tridiagonal is more efficient for tridiagonal matrices.
from scipy.linalg import eigh_tridiagonal

solve()