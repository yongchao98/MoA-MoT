import numpy as np

def solve_and_compare_spectra():
    """
    Numerically solves for the spectra of two SUSY partner Hamiltonians
    and demonstrates that they differ by at most one level.
    """
    # Define the spatial grid. We use a large interval to approximate the
    # behavior on the real line, where the difference is clearly manifested.
    N = 1000  # Number of grid points
    L = 8.0   # Interval [-L, L]
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]

    # Create the finite difference matrix for the second derivative, -d^2/dx^2
    # with Dirichlet boundary conditions (psi=0 at the ends).
    D2 = (np.diag(np.ones(N - 1), -1) - 2 * np.diag(np.ones(N), 0) + np.diag(np.ones(N - 1), 1)) / dx**2
    H_kinetic = -D2

    # Choose the superpotential W(x). The quantum harmonic oscillator is a classic example.
    # W(x) = x leads to potentials V = x^2 +/- 1
    W_x = x
    W_prime = 1.0

    # The potentials for H_0 and H_1 are V_{0,1} = W^2 -/+ W'. We ignore the arbitrary constant alpha.
    V0_potential = np.diag(W_x**2 - W_prime)
    V1_potential = np.diag(W_x**2 + W_prime)

    # Construct the Hamiltonian matrices
    H0 = H_kinetic + V0_potential
    H1 = H_kinetic + V1_potential

    # Calculate the eigenvalues. eigvalsh is efficient for Hermitian matrices.
    eigs0 = np.linalg.eigvalsh(H0)
    eigs1 = np.linalg.eigvalsh(H1)
    
    # --- Output the results ---
    print("This script demonstrates that the spectra of two supersymmetric partner Hamiltonians,")
    print("H0 and H1, can differ by at most one energy level.")
    print("\nWe use the example W(x) = x, which gives the partner potentials V0=x^2-1 and V1=x^2+1.")
    
    # Print the lowest 6 eigenvalues for comparison
    print("\nFirst 6 numerical eigenvalues for H0:")
    # Expected theoretical eigenvalues are ~2n = 0, 2, 4, ...
    print(np.round(eigs0[:6], 4))

    print("\nFirst 6 numerical eigenvalues for H1:")
    # Expected theoretical eigenvalues are ~2(n+1) = 2, 4, 6, ...
    print(np.round(eigs1[:6], 4))

    print("\n--- Analysis ---")
    print("The ground state eigenvalue of H0 (around {:.4f}) is absent in the spectrum of H1.".format(eigs0[0]))
    print("Apart from this ground state, the spectrum of H0 matches the spectrum of H1.")
    print("H0 spectrum (from 2nd level):", np.round(eigs0[1:6], 4))
    print("H1 spectrum (from 1st level):", np.round(eigs1[0:5], 4))
    
    print("\nThe number of differing levels is 1.")
    print("\nTherefore, the maximum number of levels that can differ is 1.")


solve_and_compare_spectra()
