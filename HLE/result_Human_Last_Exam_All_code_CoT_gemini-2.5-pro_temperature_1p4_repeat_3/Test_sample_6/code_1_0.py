import numpy as np

def count_kk_modes():
    """
    This function calculates the number of spin-2 Kaluza-Klein modes with
    mass-squared below a given threshold for the specified 5D gravitational theory.
    """
    
    # 1. Define problem parameters
    threshold = 14
    # The warp factor is A(x) = sin(x) + 4*cos(x).
    # The potential for the Schrödinger equation is V(x) = (3/2)A''(x) + (9/4)(A'(x))^2.
    # The Fourier coefficients V_k for V(x) can be calculated analytically.
    # V(x) = (153/8) - 6*cos(x) - (135/8)*cos(2*x) - (3/2)*sin(x) - 9*sin(2*x)
    # V_0 = 153/8 = 19.125
    # V_1 = -3 + 0.75j, V_-1 = -3 - 0.75j
    # V_2 = -135/16 + 4.5j = -8.4375 + 4.5j, V_-2 = -8.4375 - 4.5j
    # Other V_k are zero.
    V_coeffs = {
        0: 19.125,
        1: -3 + 0.75j,
        -1: -3 - 0.75j,
        2: -8.4375 + 4.5j,
        -2: -8.4375 - 4.5j
    }

    # 2. Set up the Hamiltonian matrix in Fourier space.
    # We truncate the Fourier series at a mode N. A larger N gives more accurate results.
    # N=20 is sufficient for the lowest eigenvalues to converge.
    N = 20
    dim = 2 * N + 1
    H = np.zeros((dim, dim), dtype=complex)
    
    # The basis vectors are e^{inx} for n = -N, ..., N.
    # The matrix elements are H_mn = n^2 * delta_mn + V_{m-n}.
    for m_idx in range(dim):
        for n_idx in range(dim):
            m = m_idx - N
            n = n_idx - N
            
            # Kinetic part (diagonal)
            if m == n:
                H[m_idx, n_idx] += n**2
            
            # Potential part (from V_{m-n})
            H[m_idx, n_idx] += V_coeffs.get(m - n, 0)

    # 3. Find the eigenvalues of the Hermitian matrix H.
    # np.linalg.eigvalsh is efficient for Hermitian matrices and returns sorted real eigenvalues.
    eigenvalues = np.linalg.eigvalsh(H)

    # 4. Count the eigenvalues below the threshold.
    low_eigenvalues = eigenvalues[eigenvalues < threshold]
    count = len(low_eigenvalues)

    # 5. Print the results in a clear format, including the "final equation".
    print(f"The mass-squared eigenvalues (m^2) of the spin-2 Kaluza-Klein modes are determined by solving a Schrödinger equation.")
    print(f"We are counting the number of eigenvalues below the threshold of {threshold}.")
    
    print("\nThe eigenvalues found below 14 are:")
    if count == 0:
        print("None")
    else:
        # We represent each found eigenvalue with a '1' in the sum.
        equation_str = " + ".join(["1"] * count)
        # Print the numerical values of the eigenvalues found.
        eigenvalue_values_str = ", ".join([f"{val:.4f}" for val in low_eigenvalues])
        print(eigenvalue_values_str)

        print("\nCounting each of these eigenvalues gives the total number:")
        print(f"{equation_str} = {count}")

    return count

if __name__ == "__main__":
    final_count = count_kk_modes()
    print(f"\n<<< {final_count} >>>")