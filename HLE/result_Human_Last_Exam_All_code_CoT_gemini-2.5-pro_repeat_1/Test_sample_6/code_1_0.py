import numpy as np

def solve_kk_masses():
    """
    This function calculates the Kaluza-Klein mass-squared eigenvalues for a
    5D gravitational theory and counts how many are below a given threshold.

    The problem is converted to a 1D Schrödinger equation on a circle, which is
    then solved numerically by diagonalizing the Hamiltonian in a truncated
    Fourier basis.
    """

    # The potential V(x) is derived from A(x) = sin(x) + 4cos(x).
    # V(x) = (9/4)*(A'(x))^2 + (3/2)*A''(x)
    # The Fourier coefficients v_k of V(x) = sum(v_k * exp(i*k*x)) are calculated to be:
    v0 = 153.0 / 8.0
    v1 = -3.0 + 0.75j
    v2 = -135.0 / 16.0 + 4.5j
    v = {
        0: v0,
        1: v1,
        -1: np.conj(v1),
        2: v2,
        -2: np.conj(v2)
    }

    # Set the truncation for the Fourier basis. A larger N_max gives more accurate results.
    # N_max = 30 is sufficient for the lowest eigenvalues to converge.
    N_max = 30
    matrix_size = 2 * N_max + 1

    # Construct the Hamiltonian matrix H in the Fourier basis {e^{inx}}.
    # The matrix elements are H_{m,n} = m^2 * delta_{m,n} + v_{m-n}.
    H = np.zeros((matrix_size, matrix_size), dtype=np.complex128)

    for i in range(matrix_size):
        m = i - N_max  # Corresponds to Fourier mode e^{imx}
        for j in range(matrix_size):
            n = j - N_max  # Corresponds to Fourier mode e^{inx}
            
            # The kinetic term -d^2/dx^2 becomes m^2 on the diagonal
            if m == n:
                H[i, j] += m**2
            
            # The potential term V(x) contributes v_{m-n}
            k = m - n
            if k in v:
                H[i, j] += v[k]

    # The Hamiltonian H is Hermitian, so its eigenvalues are real.
    # We use np.linalg.eigvalsh which is optimized for Hermitian matrices.
    eigenvalues = np.linalg.eigvalsh(H)

    # Filter the eigenvalues that are below the threshold of 14.
    threshold = 14
    below_threshold_eigs = eigenvalues[eigenvalues < threshold]

    print("The eigenvalues (mass-squared values) below 14 are:")
    # The final equation is the list of eigenvalues. We print each one.
    for eig in below_threshold_eigs:
        # The lowest eigenvalue corresponds to the massless 4D graviton and should be 0.
        # Numerical precision might result in a very small number instead.
        if np.abs(eig) < 1e-9:
            print(0.0)
        else:
            print(eig)
    
    print(f"\nThe number of eigenvalues below {threshold} is:")
    print(len(below_threshold_eigs))

solve_kk_masses()