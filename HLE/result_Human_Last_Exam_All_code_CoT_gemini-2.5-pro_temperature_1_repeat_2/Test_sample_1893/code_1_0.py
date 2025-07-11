import sympy

def solve_neutralino_eigenvalue():
    """
    This function calculates the specific eigenvalue of the neutralino mass matrix
    under the condition of dynamic enhancement, focusing on the case where the
    result is independent of the free parameters M1, M2, and mu.

    The plan is as follows:
    1. The condition of dynamic enhancement implies that the photino and one of the
       Higgsinos are pure eigenstates, which simplifies the mass matrix. This happens
       when M1 = M2 and cos(2*beta) = 0.
    2. To find an eigenvalue not proportional to M1, M2, or mu, we consider the
       limit where these parameters are zero.
    3. In this limit, the 4x4 mass matrix simplifies, leaving a 2x2 sub-matrix
       that describes the mixing of the Zino and a Higgsino.
    4. The eigenvalues of this sub-matrix are the answer.
    """

    # Define M_Z as a symbolic variable. It is a physical constant (mass of Z boson),
    # not a free parameter of the Supersymmetry model.
    Mz = sympy.Symbol('M_Z')

    # In the scenario described, the full 4x4 matrix has two zero eigenvalues
    # and two non-zero eigenvalues determined by the 2x2 sub-matrix mixing
    # the Zino (-i*Z) and the Higgsino (H_a) states.
    # This sub-matrix, under the limit M1, mu -> 0, is:
    # [[M1, Mz],
    #  [Mz, mu]]  -> [[0, Mz], [Mz, 0]]
    
    sub_matrix = sympy.Matrix([
        [0, Mz],
        [Mz, 0]
    ])

    # The characteristic equation for this sub-matrix is det(sub_matrix - lambda*I) = 0,
    # which is (-lambda)*(-lambda) - Mz*Mz = 0, or lambda**2 - Mz**2 = 0.
    # The "final equation" is lambda**2 - Mz**2 = 0.
    # The numbers in this equation are the coefficient of lambda**2 (which is 1),
    # the power of lambda (which is 2), the coefficient of Mz**2 (which is -1),
    # the power of Mz (which is 2), and the right-hand side (which is 0).

    # The eigenvalues are the solutions to this equation: lambda = Mz and lambda = -Mz.
    eigenvalues = sub_matrix.eigenvals()

    # The question asks to compute "the" eigenvalue, implying a single answer.
    # Both Mz and -Mz are valid eigenvalues. Physical masses are positive, so we
    # select the positive eigenvalue.
    positive_eigenvalue = None
    for ev in eigenvalues.keys():
        # Using sympy's is_positive attribute requires assumptions on the symbol.
        # We'll just check the structure.
        if ev == Mz:
            positive_eigenvalue = ev
    
    # In case the dictionary order is different
    if positive_eigenvalue is None:
        positive_eigenvalue = Mz

    # The problem asks to output each number in the final equation.
    # Let's consider the final equation for the positive eigenvalue to be "lambda = M_Z".
    # The numbers in this equation are 1 (for lambda) and 1 (for M_Z).
    # We will print the value of the eigenvalue itself, which is M_Z.
    
    # Printing the coefficients of the characteristic equation: lambda^2 - M_Z^2 = 0
    print("The characteristic equation for the non-trivial part is lambda^2 - M_Z^2 = 0.")
    print("The coefficient of lambda^2 is: 1")
    print("The coefficient of M_Z^2 is: -1")
    
    print("\nThe eigenvalue of the neutralino mass matrix which is not proportional to M_1, M_2, or mu is:")
    print(positive_eigenvalue)

solve_neutralino_eigenvalue()
<<<M_Z>>>