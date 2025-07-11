import numpy as np

def count_kk_modes():
    """
    This function calculates the number of spin-2 KK mode eigenvalues below a threshold
    for a given warped compactification.
    """
    # The problem is to find the number of eigenvalues of a Schrödinger operator
    # H = -d^2/dx^2 + V(x) that are less than 14.
    # The potential is V(x) = (9/4)*(A'(x))^2 + (3/2)*A''(x)
    # with A(x) = sin(x) + 4*cos(x).
    # The problem is on a circle x in [0, 2*pi], so we use a Fourier basis.

    # 1. Define a cutoff for the Fourier basis {e^(i*n*x)}
    n_max = 50
    # The modes are indexed from -n_max to n_max
    n_range = np.arange(-n_max, n_max + 1)
    # Total number of basis states (size of the matrix)
    N = len(n_range)

    # 2. Calculate the Fourier coefficients of the potential V(x).
    # V(x) = 153/8 - 6*cos(x) - 3/2*sin(x) - 135/8*cos(2x) - 18*sin(2x)
    # The coefficients V_k are such that V(x) = sum_k V_k * exp(i*k*x).
    # Note that V_k corresponds to the coefficient of exp(i*(-k)*x) in the standard
    # physics convention for Fourier series, which is what we need for H_mn = V_{m-n}.
    V0 = 153.0 / 8.0
    V1 = -3.0 - 0.75j        # H_{m, m-1} element
    V_m1 = -3.0 + 0.75j      # H_{m, m+1} element
    V2 = -135.0 / 16.0 - 9.0j  # H_{m, m-2} element
    V_m2 = -135.0 / 16.0 + 9.0j  # H_{m, m+2} element

    # 3. Construct the Hamiltonian matrix H_mn = n^2 * delta_mn + V_{m-n}
    # H is a pentadiagonal matrix.
    
    # The main diagonal contains the kinetic term n^2 and the constant part of the potential V0.
    diag0 = n_range**2 + V0
    
    # Construct the matrix from its diagonals.
    # We start with the main diagonal.
    H = np.diag(diag0)
    # Add the other four non-zero diagonals.
    H = H + np.diag(np.full(N - 1, V1, dtype=np.complex128), k=-1)
    H = H + np.diag(np.full(N - 1, V_m1, dtype=np.complex128), k=1)
    H = H + np.diag(np.full(N - 2, V2, dtype=np.complex128), k=-2)
    H = H + np.diag(np.full(N - 2, V_m2, dtype=np.complex128), k=2)

    # 4. Find the eigenvalues of the Hermitian matrix H.
    # numpy.linalg.eigh is specialized for Hermitian matrices and returns real eigenvalues.
    eigenvalues = np.linalg.eigh(H)[0]

    # 5. Count the number of eigenvalues below the specified threshold.
    threshold = 14
    count = np.sum(eigenvalues < threshold)

    # The problem asks for the final answer to be printed.
    # The final "equation" is simply the result of this counting.
    # The numbers involved are the threshold and the count.
    print(f"Number of eigenvalues below {threshold}: {count}")
    # The final equation is essentially: count = 4
    # The numbers in this equation are 4. Let's just print the count.
    print("Final answer:")
    print(int(count))

count_kk_modes()