import numpy as np
from scipy.linalg import eigh_tridiagonal

def solve_schrodinger(potential_func, x_grid):
    """
    Solves the 1D time-independent Schrödinger equation using finite differences.
    Returns the first few eigenvalues.
    """
    N = len(x_grid)
    dx = x_grid[1] - x_grid[0]
    
    # Potential energy on the diagonal
    diagonal = potential_func(x_grid) + 2.0 / (dx**2)
    
    # Kinetic energy on the off-diagonal
    off_diagonal = -1.0 / (dx**2) * np.ones(N - 1)
    
    # Solve the tridiagonal eigenvalue problem
    # eigh_tridiagonal is much faster for this specific type of matrix.
    eigenvalues, _ = eigh_tridiagonal(diagonal, off_diagonal)
    
    return eigenvalues

def main():
    """
    Main function to define SUSY partner Hamiltonians and compare their spectra.
    """
    # --- 1. Define the SUSY model ---
    # We choose a superpotential W(x) and an energy shift alpha.
    # The example is based on the Quantum Harmonic Oscillator.
    # W(x) = -x  (Note: dW/dx = -1)
    # alpha = -1
    def W(x):
        return -x
        
    def W_prime(x):
        return -1.0 * np.ones_like(x)

    alpha = -1.0

    # --- 2. Define the partner potentials V0 and V1 ---
    # According to the problem statement and the factorization:
    # H_0 = -d^2/dx^2 + V_0(x) where V_0(x) = W(x)^2 - W'(x) - alpha
    # H_1 = -d^2/dx^2 + V_1(x) where V_1(x) = W(x)^2 + W'(x) - alpha
    def V0(x):
        return W(x)**2 - W_prime(x) - alpha

    def V1(x):
        return W(x)**2 + W_prime(x) - alpha

    # This choice leads to:
    # V0(x) = (-x)^2 - (-1) - (-1) = x^2 + 2
    # V1(x) = (-x)^2 + (-1) - (-1) = x^2
    # So we compare H_0 = -d^2/dx^2 + x^2 + 2 and H_1 = -d^2/dx^2 + x^2

    # --- 3. Numerical parameters ---
    N = 2001  # Number of grid points
    L = 8.0   # Grid extends from -L to L
    x_grid = np.linspace(-L, L, N)
    
    # --- 4. Solve for eigenvalues ---
    print("Calculating eigenvalues for H_0 and H_1...")
    eigenvalues_H0 = solve_schrodinger(V0, x_grid)
    eigenvalues_H1 = solve_schrodinger(V1, x_grid)

    # --- 5. Compare spectra ---
    num_levels_to_show = 5
    print("\nFirst {} numerical eigenvalues:".format(num_levels_to_show))
    print("--------------------------------")
    print("H_0 (V=x^2+2):   H_1 (V=x^2):")
    for i in range(num_levels_to_show):
        print(f"E{i}: {eigenvalues_H0[i]:<14.4f} E{i}: {eigenvalues_H1[i]:<14.4f}")
    
    # Analytically, for H_QHO = -d^2/dx^2 + x^2, eigenvalues are 1, 3, 5, 7, ...
    # So, spec(H_1) = {1, 3, 5, 7, ...}
    # and spec(H_0) = {1+2, 3+2, 5+2, 7+2, ...} = {3, 5, 7, 9, ...}
    
    print("\n--- Analysis ---")
    print("The theoretical spectrum for H_1 is {1, 3, 5, 7, ...}")
    print("The theoretical spectrum for H_0 is {3, 5, 7, 9, ...}")
    print("\nThe spectrum of H_0 is almost identical to the spectrum of H_1,")
    print("but H_1 has one energy level (E=1) that H_0 does not possess.")
    print("This level E=1 corresponds to -alpha.")
    print(f"Energy shift alpha = {alpha}, so the special energy level is E = -alpha = {-alpha}")

    # The maximum number of levels that can differ is determined by theory.
    # The code illustrates the case where this maximum is achieved.
    max_diff = 1
    
    print("\n--- Conclusion ---")
    # This prints the required final answer.
    print(f"The final equation H_0 = L^+ L - ({alpha}) and H_1 = L L^+ - ({alpha}) shows that the spectra can differ at most at the energy E = {-alpha}.")
    print(f"The maximum number of levels of the spectrum that can differ is: {max_diff}")

if __name__ == '__main__':
    main()
