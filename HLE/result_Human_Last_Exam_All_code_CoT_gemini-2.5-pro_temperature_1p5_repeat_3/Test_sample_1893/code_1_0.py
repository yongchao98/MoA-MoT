import numpy as np

def solve_neutralino_mass():
    """
    Computes the specific eigenvalue of the neutralino mass matrix based on the problem's conditions.

    The problem describes a "dynamic enhancement" scenario where the photino and one of the Higgsinos
    decouple from the zino and the other Higgsino. This requires M1 = M2 and cos(2*beta) = 0.
    
    To find an eigenvalue not proportional to the adjustable parameters M1, M2, or mu, we
    consider the special case where M1 = M2 = 0 and mu = 0. This is a valid choice that satisfies
    the decoupling conditions and reveals eigenvalues that are physical constants.
    """
    
    # Physical constants
    M_Z = 91.1876  # Z boson mass in GeV
    sin2_theta_W = 0.2223  # Square of the sine of the Weinberg angle
    s = np.sqrt(sin2_theta_W)
    c = np.sqrt(1 - sin2_theta_W)

    # Parameters from the problem's conditions
    M1 = 0.0
    M2 = 0.0  # Condition M1 = M2
    mu = 0.0
    # Condition cos(2*beta) = 0 implies beta = pi/4 for 0 <= beta <= pi/2
    beta = np.pi / 4.0 

    # Construct the neutralino mass matrix M_N
    M_N = np.zeros((4, 4))
    
    M_N[0, 0] = M1 * c**2 + M2 * s**2
    M_N[0, 1] = (M2 - M1) * s * c
    M_N[1, 0] = (M2 - M1) * s * c
    M_N[1, 1] = M1 * s**2 + M2 * c**2
    M_N[1, 2] = M_Z
    M_N[2, 1] = M_Z
    M_N[2, 2] = mu * np.sin(2 * beta)
    M_N[2, 3] = -mu * np.cos(2 * beta)
    M_N[3, 2] = -mu * np.cos(2 * beta)
    M_N[3, 3] = -mu * np.sin(2 * beta)

    print("Under the specified conditions (M1=M2=0, mu=0, beta=pi/4), the neutralino mass matrix becomes:")
    # The command below prints all the numbers in the matrix
    print(np.round(M_N, 4))

    # Calculate the eigenvalues of the matrix
    eigenvalues = np.linalg.eigvals(M_N)

    print("\nThe four eigenvalues of this matrix are:")
    print(np.round(eigenvalues, 4))

    # Find the positive, non-zero eigenvalue
    result = 0.0
    for eig in eigenvalues:
        if eig > 1e-6: # Check for a positive value, allowing for float inaccuracy
            result = eig
            break
            
    print(f"\nThe eigenvalue not proportional to M1, M2, or mu is the mass of the Z boson:")
    print(f"{result:.4f}")

solve_neutralino_mass()
<<<91.1876>>>