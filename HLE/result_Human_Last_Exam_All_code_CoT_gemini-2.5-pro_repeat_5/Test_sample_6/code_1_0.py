import numpy as np

def solve_kk_masses():
    """
    This function calculates the Kaluza-Klein mass spectrum for a given
    warp factor and counts the number of modes below a specified threshold.
    """
    # 1. Set up the numerical grid
    N = 1000  # Number of grid points for high accuracy
    x, dx = np.linspace(0, 2 * np.pi, N, endpoint=False, retstep=True)

    # 2. Define the warp factor and calculate the potential V(x)
    A = np.sin(x) + 4 * np.cos(x)
    A_prime = np.cos(x) - 4 * np.sin(x)
    # The potential is V(x) = (3/2)A''(x) + (9/4)A'(x)^2, where A''(x) = -A(x)
    V = - (3.0 / 2.0) * A + (9.0 / 4.0) * (A_prime**2)

    # 3. Construct the Hamiltonian matrix for the Schrödinger operator
    # H = -d^2/dx^2 + V(x) using finite differences with periodic BCs.
    H = np.zeros((N, N))
    
    # Main diagonal
    main_diag = 2 / dx**2 + V
    np.fill_diagonal(H, main_diag)
    
    # Off-diagonals
    off_diag = -1 / dx**2
    np.fill_diagonal(H[1:], off_diag)
    np.fill_diagonal(H[:, 1:], off_diag)
    
    # Corner elements for periodic boundary conditions
    H[0, N - 1] = off_diag
    H[N - 1, 0] = off_diag

    # 4. Calculate the eigenvalues (m^2) of the Hamiltonian matrix.
    # We use eigvalsh because H is real and symmetric.
    eigenvalues_m_squared = np.linalg.eigvalsh(H)

    # 5. Count the number of modes with mass m < 14, i.e., m^2 < 14^2
    mass_threshold = 14.0
    mass_squared_threshold = mass_threshold**2
    
    modes_below_threshold = eigenvalues_m_squared[eigenvalues_m_squared < mass_squared_threshold]
    
    count = len(modes_below_threshold)

    # 6. Output the results as requested
    print(f"The eigenvalues (m^2) below the threshold of {mass_squared_threshold:.2f} are:")
    # The "final equation" is the sum of the modes found. We print each one.
    for val in modes_below_threshold:
        print(f"{val:.4f}")
        
    print(f"\nThe final count of eigenvalues below the threshold is: {count}")

solve_kk_masses()
<<<27>>>