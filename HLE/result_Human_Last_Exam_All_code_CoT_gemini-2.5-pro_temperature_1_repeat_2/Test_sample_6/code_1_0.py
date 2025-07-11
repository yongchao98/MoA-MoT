import numpy as np

def count_kk_modes():
    """
    Calculates the number of spin-2 KK mode eigenvalues below 14.

    The problem is mapped to a Schrödinger equation -phi'' + V(x)phi = m^2 phi
    with a potential derived from the warp factor A(x).
    V(x) = 4*(A'(x))^2 + 2*A''(x), with A(x) = sin(x) + 4cos(x).
    This gives V(x) = 34 - 8*cos(x) - 2*sin(x) - 30*cos(2*x) - 16*sin(2*x).

    We solve this using a Fourier basis expansion, which turns the problem
    into a matrix eigenvalue problem Hc = m^2 c. The matrix elements are
    H_{mn} = n^2 * delta_{mn} + V_hat_{m-n}, where V_hat are the Fourier
    coefficients of the potential.
    """

    # Set the truncation for the Fourier basis. A larger N_max improves accuracy.
    # N_max = 30 is sufficient for the low-lying eigenvalues to converge.
    N_max = 30

    # The Fourier coefficients V_hat_k of the potential V(x) are derived from its
    # expression in terms of complex exponentials.
    # V(x) = 34 + (-4+i)e^{-ix} + (-4-i)e^{ix} + (-15+8i)e^{-i2x} + (-15-8i)e^{i2x}
    # These numbers define the off-diagonal elements of the Hamiltonian matrix.
    V_hat = {
        0: 34.0,
        1: -4.0 - 1.0j,
        -1: -4.0 + 1.0j,
        2: -15.0 - 8.0j,
        -2: -15.0 + 8.0j
    }

    # The size of the Hamiltonian matrix is (2*N_max + 1) x (2*N_max + 1)
    dim = 2 * N_max + 1
    H = np.zeros((dim, dim), dtype=np.complex128)

    # Populate the Hamiltonian matrix H_{mn}
    # The matrix indices i, j correspond to Fourier modes m, n from -N_max to N_max
    for i in range(dim):
        m = i - N_max
        for j in range(dim):
            n = j - N_max
            
            # The kinetic term n^2 is on the diagonal
            if m == n:
                H[i, j] += n**2
            
            # The potential term V_hat_{m-n} provides the off-diagonal elements
            k = m - n
            if k in V_hat:
                H[i, j] += V_hat[k]

    # The Hamiltonian matrix is Hermitian, so its eigenvalues (m^2) are real.
    # np.linalg.eigvalsh is efficient for Hermitian matrices and returns sorted eigenvalues.
    eigenvalues = np.linalg.eigvalsh(H)

    # Count how many eigenvalues are below the threshold of 14.
    count = np.sum(eigenvalues < 14)

    print(count)

# Execute the function to find the answer.
count_kk_modes()