import numpy as np

def count_kk_modes():
    """
    This script calculates the number of spin-2 Kaluza-Klein modes with
    mass-squared below 14 for the given 5D warped compactification.

    The problem is equivalent to finding the number of eigenvalues below 14
    for the Schrödinger operator H = -d^2/dx^2 + V(x) on a circle.

    The potential V(x) is derived from the warp factor A(x) = sin(x) + 4*cos(x)
    and is found to be:
    V(x) = 34 - 8*cos(x) - 2*sin(x) - 30*cos(2x) - 16*sin(2x)
    
    We solve this by diagonalizing the Hamiltonian in a truncated Fourier basis.
    """

    # Set the cutoff for the Fourier modes. A value of 50 is sufficient for convergence.
    N = 50
    size = 2 * N + 1

    # An array of the integer mode numbers n, from -N to N.
    n_vals = np.arange(-N, N + 1)

    # The Fourier coefficients V_k of the potential V(x) are derived from its
    # trigonometric expression. V_k is the coefficient of exp(i*k*x).
    # For V(x) = 34 - 8*cos(x) - 2*sin(x) - 30*cos(2x) - 16*sin(2x):
    # V_0 = 34
    # V_1 corresponds to (-8*cos(x) - 2*sin(x)) -> -4+j
    # V_2 corresponds to (-30*cos(2x) - 16*sin(2x)) -> -15+8j
    # and V_k = V_{-k}^*
    V_coeffs = {
        0: 34.0,           # Constant term
        1: -4.0 + 1.0j,    # Coefficient for exp(i*x)
        -1: -4.0 - 1.0j,   # Coefficient for exp(-i*x)
        2: -15.0 + 8.0j,   # Coefficient for exp(2*i*x)
        -2: -15.0 - 8.0j    # Coefficient for exp(-2*i*x)
    }

    # Construct the Hamiltonian matrix H_mn = n^2 * delta_mn + V_{m-n}
    H = np.zeros((size, size), dtype=complex)

    # Populate the matrix
    for i, m in enumerate(n_vals):
        for j, n in enumerate(n_vals):
            # The kinetic energy part is diagonal: n^2
            if i == j:
                H[i, j] += n**2
            # The potential energy part depends on the difference of modes m-n
            k = m - n
            if k in V_coeffs:
                H[i, j] += V_coeffs[k]

    # The Hamiltonian H is Hermitian. We use np.linalg.eigh, which is optimized
    # for Hermitian matrices and returns real eigenvalues.
    eigenvalues = np.linalg.eigh(H)[0]

    # The threshold for the eigenvalues (squared masses) is 14.
    threshold = 14

    # Count how many eigenvalues are below the threshold.
    count = np.sum(eigenvalues < threshold)

    print(f"The number of eigenvalues below {threshold} is: {count}")

# Execute the function to find the answer.
count_kk_modes()