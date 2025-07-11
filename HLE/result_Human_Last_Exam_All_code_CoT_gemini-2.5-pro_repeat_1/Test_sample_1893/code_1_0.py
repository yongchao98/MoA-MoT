import numpy as np

def solve_neutralino_eigenvalue():
    """
    Calculates the specified neutralino mass matrix eigenvalue.

    The problem asks for the eigenvalue that is not proportional to M_1, M_2, or mu
    under the condition of dynamic enhancement. This occurs in the specific scenario where
    M_1 = M_2 = mu = 0. In this case, the mass matrix simplifies, and the
    characteristic equation for the eigenvalues (lambda) becomes:
    lambda^2 * (lambda^2 - M_Z^2) = 0.
    The eigenvalues are 0, 0, M_Z, and -M_Z. The non-zero eigenvalues are
    determined by the Z boson mass, M_Z. We will return the positive value.
    """

    # The mass of the Z boson in GeV. This is a well-measured physical constant.
    M_Z = 91.1876

    # Under the condition M_1 = M_2 = mu = 0, the neutralino mass matrix becomes very sparse.
    # In the given basis ( -i*gamma, -i*Z, H_a, H_b ), the matrix is:
    # M_N = [[0, 0,  0,  0],
    #        [0, 0, M_Z, 0],
    #        [0, M_Z, 0,  0],
    #        [0, 0,  0,  0]]
    M_N = np.array([
        [0., 0., 0., 0.],
        [0., 0., M_Z, 0.],
        [0., M_Z, 0., 0.],
        [0., 0., 0., 0.]
    ])

    # We can find the eigenvalues of this matrix.
    eigenvalues = np.linalg.eigvals(M_N)

    # The eigenvalue that fits the description is the positive, non-zero one.
    result_eigenvalue = 0.0
    for val in eigenvalues:
        if val > 1e-9:  # Use tolerance for float comparison
            result_eigenvalue = val
            break
            
    # The final equation for the non-zero eigenvalues is lambda^2 = M_Z^2.
    # The problem asks to output each number in the final equation. Here, this is M_Z.
    print(f"The value of M_Z in the equation is: {M_Z}")
    print(f"The computed positive eigenvalue is: {result_eigenvalue}")

solve_neutralino_eigenvalue()