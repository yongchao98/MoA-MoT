import sympy as sp

def solve_eigenvalue():
    """
    This function computes the specific eigenvalue of the neutralino mass matrix
    under the conditions described in the problem.
    """
    # Define symbolic variables for the parameters in the matrix.
    # M1, M2, mu are adjustable parameters. MZ is the Z boson mass (a constant).
    # theta_W is the Weinberg angle. beta is a mixing angle.
    M1, M2, mu, beta, MZ = sp.symbols('M1 M2 mu beta MZ', real=True)
    sw, cw = sp.symbols('sw cw', real=True, positive=True)

    # The problem describes a "dynamic enhancement" scenario. This occurs when
    # the subspace spanned by gamma-tilde (basis vector 1) and Hb-tilde (basis vector 4)
    # is invariant under the mass matrix transformation.
    # This implies that when the matrix acts on a vector (a, 0, 0, d), the result
    # must be of the form (a', 0, 0, d').
    # Let's inspect the matrix multiplication:
    # M * [a, 0, 0, d]^T = [..., a*(M2 - M1)*sw*cw, d*(-mu*cos(2*beta)), ...]^T
    # For the 2nd and 3rd components to be zero for any 'a' and 'd', their coefficients
    # must be zero.
    # Condition 1: (M2 - M1) * sw * cw = 0. Since sw and cw are non-zero, this means M1 = M2.
    # Condition 2: -mu * cos(2*beta) = 0. This means either mu = 0 or cos(2*beta) = 0.

    # We are looking for an eigenvalue that is NOT proportional to M1, M2, or mu.
    # This suggests the eigenvalue is a constant. Let's analyze the two possibilities
    # from Condition 2.
    #
    # Case 1: cos(2*beta) = 0 and mu != 0. The eigenvalues in this case will depend on M1 and mu.
    # There is no constant eigenvalue in this case.
    #
    # Case 2: mu = 0. Let's see if this gives a constant eigenvalue.
    # If we set mu = 0 and M1 = M2, the decoupling conditions are met.

    # Let's construct the matrix with M1=M2 and mu=0.
    # Note that sw**2 + cw**2 = 1.
    M_N_special = sp.Matrix([
        [M1, 0, 0, 0],
        [0, M1, MZ, 0],
        [0, MZ, 0, 0],
        [0, 0, 0, 0]
    ])

    # The eigenvalues of this matrix can be found easily. Due to the block-diagonal
    # structure and the row of zeros, we can see the eigenvalues by inspection.
    # The characteristic polynomial is det(M_N_special - lambda*I).
    # (M1 - lambda) * (-lambda) * det([[M1-lambda, MZ], [MZ, -lambda]]) = 0
    # The solutions are lambda=M1, lambda=0, and the roots of (lambda^2 - M1*lambda - MZ^2 = 0).
    
    # The eigenvalues are:
    # 0
    # M1
    # (M1 + sqrt(M1**2 + 4*MZ**2))/2
    # (M1 - sqrt(M1**2 + 4*MZ**2))/2
    
    # Of these four eigenvalues, only one is a constant that does not depend on M1, M2, or mu.
    constant_eigenvalue = 0

    print("The scenario described requires M1=M2 and mu=0 to produce an eigenvalue independent of the adjustable parameters.")
    print("Under these conditions, one of the eigenvalues of the neutralino mass matrix is 0.")
    print("The final equation for the requested eigenvalue is: Eigenvalue = 0.")
    print("The computed numerical value of the eigenvalue is:")
    print(constant_eigenvalue)

solve_eigenvalue()