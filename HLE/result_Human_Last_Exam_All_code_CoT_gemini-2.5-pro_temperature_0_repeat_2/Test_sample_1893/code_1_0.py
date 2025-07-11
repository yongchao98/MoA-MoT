import numpy as np

def solve_neutralino_mass():
    """
    This function calculates the specific eigenvalue of the neutralino mass matrix
    as requested by the user.
    """
    # The mass of the Z boson in GeV. This is a physical constant.
    M_Z = 91.1876

    # The problem asks for an eigenvalue under several conditions.
    # 1. Dynamic enhancement conditions: M1 = M2 and beta = pi/4.
    # 2. The eigenvalue should not be proportional to M1, M2, or mu.
    # We interpret this as finding the eigenvalue in the limit where the
    # adjustable mass parameters M1 and mu are zero.
    M1 = 0.0
    mu = 0.0

    # Under these conditions, the neutralino mass matrix simplifies to:
    # M = [[M1, 0,   0,   0],
    #      [0,  M1,  M_Z, 0],
    #      [0,  M_Z, mu,  0],
    #      [0,  0,   0,  -mu]]
    #
    # Substituting M1 = 0 and mu = 0, we get the final matrix.
    M_final = np.array([
        [0.0, 0.0, M_Z, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [M_Z, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0]
    ])
    
    # The provided matrix is in the basis (g, Z, Ha, Hb).
    # Let's reorder it to match the problem's basis (-i*g, -i*Z, Ha, Hb)
    # The matrix from the problem with M1=0, mu=0, beta=pi/4 is:
    M_final_problem_basis = np.array([
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, M_Z, 0.0],
        [0.0, M_Z, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0]
    ])


    # The characteristic equation for this matrix is det(M - lambda*I) = 0,
    # which simplifies to: lambda^2 * (lambda^2 - M_Z^2) = 0.
    # The numbers in this equation are 1 (for lambda^4), -M_Z^2, and 0.
    
    # We can find the eigenvalues numerically.
    eigenvalues = np.linalg.eigvals(M_final_problem_basis)

    # We are looking for the eigenvalue that is not zero.
    # The eigenvalues will be 0, 0, M_Z, and -M_Z.
    # We select the positive mass eigenvalue.
    result_eigenvalue = 0.0
    for val in eigenvalues:
        # Use a small tolerance for floating point comparison to find the positive value
        if val > 1e-9:
            result_eigenvalue = val
            break
            
    print("The characteristic equation is: lambda^4 - (M_Z^2 * lambda^2) = 0")
    print(f"The numbers in this equation are 1, -{M_Z**2:.4f}, and 0.")
    print(f"The positive, non-zero eigenvalue is the mass of the Z boson, M_Z.")
    print("\nThe computed eigenvalue is:")
    print(result_eigenvalue)

solve_neutralino_mass()