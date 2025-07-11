import numpy as np

def solve_kk_masses():
    """
    Calculates the number of Kaluza-Klein spin-2 eigenvalues below a threshold.
    """
    # Set up the parameters for the calculation.
    # N_max determines the truncation size of the Fourier basis {-N_max, ..., N_max}.
    # A larger value gives more accurate results for the low-lying eigenvalues.
    N_max = 30
    matrix_size = 2 * N_max + 1
    threshold = 14

    # The potential V(x) is derived from the warp factor A(x) = sin(x) + 4cos(x).
    # V(x) = (9/4)*(A')^2 + (3/2)*A''.
    # We compute the Fourier coefficients V_m of V(x) = sum(V_m * exp(i*m*x)).
    # For the given A(x), only V0, V1, V-1, V2, V-2 are non-zero.
    V0 = -153.0 / 8.0
    V1 = -3.0 - 0.75j
    V_minus_1 = np.conj(V1)
    V2 = -135.0 / 16.0 - 4.5j
    V_minus_2 = np.conj(V2)
    
    # Construct the Hamiltonian matrix H in the Fourier basis.
    # The matrix element H_nk = <n|H|k> = k^2 * delta_nk + V_{n-k}.
    # The matrix is indexed from 0 to 2*N_max, corresponding to k from -N_max to N_max.
    H = np.zeros((matrix_size, matrix_size), dtype=complex)

    for row_idx in range(matrix_size):
        # Set diagonal element: k=n, so V_{n-n} = V0
        k = row_idx - N_max
        H[row_idx, row_idx] = k**2 + V0

        # Set off-diagonal elements
        # H[n,k] has V_{n-k}
        # H[row_idx, col_idx] has V_{row_idx - col_idx}
        
        # V1 for n-k=1 => col_idx = row_idx - 1
        if row_idx - 1 >= 0:
            H[row_idx, row_idx - 1] = V1
            
        # V-1 for n-k=-1 => col_idx = row_idx + 1
        if row_idx + 1 < matrix_size:
            H[row_idx, row_idx + 1] = V_minus_1

        # V2 for n-k=2 => col_idx = row_idx - 2
        if row_idx - 2 >= 0:
            H[row_idx, row_idx - 2] = V2
            
        # V-2 for n-k=-2 => col_idx = row_idx + 2
        if row_idx + 2 < matrix_size:
            H[row_idx, row_idx + 2] = V_minus_2
            
    # Diagonalize the Hermitian matrix H to find its real eigenvalues.
    # np.linalg.eigvalsh is efficient for Hermitian matrices and returns sorted eigenvalues.
    eigenvalues = np.linalg.eigvalsh(H)

    # Filter for eigenvalues below the threshold and count them.
    below_threshold_eigs = eigenvalues[eigenvalues < threshold]
    count = len(below_threshold_eigs)

    print("The eigenvalues below 14 are:")
    for ev in below_threshold_eigs:
        # We print each number that contributes to the final count.
        print(f"{ev:.4f}")
    
    print(f"\nThe total number of eigenvalues below 14 is: {count}")

solve_kk_masses()