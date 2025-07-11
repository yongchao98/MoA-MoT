import numpy as np

def solve_kk_masses():
    """
    This function calculates the Kaluza-Klein mass-squared eigenvalues for spin-2 modes
    and counts how many are below a specified threshold.
    
    The problem is mapped to a Schrödinger equation -z'' + V(x)z = m^2 z on a circle,
    which is then solved by diagonalizing the Hamiltonian in a Fourier basis.
    """
    
    # 1. Set up the problem parameters
    threshold = 14
    # Set the cutoff for the Fourier basis. A larger value gives more accurate results.
    N_max = 30
    
    # 2. Define the Fourier coefficients of the potential V(x).
    # V(x) = (9/4)*(A')^2 + (3/2)*A''
    # A(x) = sin(x) + 4*cos(x)
    # The non-zero Fourier coefficients V_m for V(x) = sum(V_m * exp(i*m*x)) are calculated analytically.
    V_coeffs = {
        0: 153.0 / 8.0,
        1: -3.0 + 0.75j,
        -1: -3.0 - 0.75j,
        2: -135.0 / 16.0 + 4.5j,
        -2: -135.0 / 16.0 - 4.5j
    }

    # 3. Construct the Hamiltonian matrix in the Fourier basis {e^{inx}}
    size = 2 * N_max + 1
    H = np.zeros((size, size), dtype=np.complex128)
    
    # The matrix elements are H_{mn} = <m|H|n> = n^2 * delta_mn + V_{m-n}
    # where m, n are momentum integers from -N_max to N_max.
    for i in range(size):
        for j in range(size):
            m = i - N_max
            n = j - N_max
            
            # Kinetic part (diagonal)
            if i == j:
                H[i, j] += n**2
            
            # Potential part
            diff = m - n
            if diff in V_coeffs:
                H[i, j] += V_coeffs[diff]

    # 4. Diagonalize the Hermitian matrix to find the eigenvalues (m^2)
    # np.linalg.eigvalsh is efficient for Hermitian matrices and returns real eigenvalues.
    eigenvalues = np.linalg.eigvalsh(H)
    
    # 5. Filter and count the eigenvalues below the threshold
    eigenvalues_below_threshold = [val for val in eigenvalues if val < threshold]
    count = len(eigenvalues_below_threshold)

    print("The mass-squared eigenvalues below 14 are:")
    for val in eigenvalues_below_threshold:
        # The lowest eigenvalue should be exactly 0, we'll format it to show that.
        if abs(val) < 1e-9:
             print(f"{0.0:.5f}")
        else:
             print(f"{val:.5f}")
    
    print(f"\nThe total number of eigenvalues below 14 is: {count}")

solve_kk_masses()
<<<5>>>