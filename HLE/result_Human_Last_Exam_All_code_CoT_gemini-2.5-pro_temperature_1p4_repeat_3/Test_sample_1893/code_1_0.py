import sympy

def solve_neutralino_mass():
    """
    Computes the special eigenvalue of the neutralino mass matrix under
    the conditions of dynamic enhancement for radiative decay.
    """
    # Define the symbolic parameters of the neutralino mass matrix
    M1, M2, mu, beta, theta_W, M_Z = sympy.symbols('M1 M2 mu beta theta_W M_Z', real=True)
    cW = sympy.cos(theta_W)
    sW = sympy.sin(theta_W)

    # Define the neutralino mass matrix in the specified basis
    # Basis: (-i*gamma_tilde, -i*Z_tilde, H_a_tilde, H_b_tilde)
    M_N = sympy.Matrix([
        [M1*cW**2 + M2*sW**2, (M2 - M1)*sW*cW, 0, 0],
        [(M2 - M1)*sW*cW, M1*sW**2 + M2*cW**2, M_Z, 0],
        [0, M_Z, mu*sympy.sin(2*beta), -mu*sympy.cos(2*beta)],
        [0, 0, -mu*sympy.cos(2*beta), -mu*sympy.sin(2*beta)]
    ])

    # Apply the conditions derived from the physics scenario:
    # 1. "pure gamma_tilde state" implies no mixing with Z_tilde -> M12 = 0 -> M1 = M2
    # 2. "dynamic enhancement" implies the one-loop coupling is dominant. This can occur
    #    if the tree-level Higgsino mass term is zero -> mu = 0
    substitutions = {M2: M1, mu: 0}
    M_N_simplified = M_N.subs(substitutions)

    # The simplified matrix becomes:
    # [[M1, 0,    0, 0],
    #  [0,  M1, M_Z, 0],
    #  [0, M_Z,   0, 0],
    #  [0,  0,    0, 0]]
    
    # Calculate the eigenvalues of the simplified matrix
    # The characteristic polynomial is det(M_N_simplified - lambda*I) = 0
    # From the last row, it's clear that one eigenvalue is 0.
    # The remaining 3x3 block has a determinant of:
    # (M1 - lambda) * ( (M1 - lambda)*(-lambda) - M_Z**2 ) = 0
    # So, one eigenvalue is M1.
    # The other two are from -lambda*(M1 - lambda) - M_Z**2 = 0
    # lambda**2 - M1*lambda + M_Z**2 = 0
    # The four eigenvalues are 0, M1, and the two roots of the quadratic equation.
    
    eigenvalues = M_N_simplified.eigenvals() # Returns a dict of {eigenvalue: multiplicity}
    
    # The question asks for the eigenvalue not proportional to M1, M2, or mu.
    # Let's find this eigenvalue from our calculated set.
    target_eigenvalue = None
    for val in eigenvalues.keys():
        # An expression is independent of a set of symbols if they don't appear in it
        if not val.has(M1, M2, mu):
            target_eigenvalue = val
            break

    print("The neutralino mass matrix under the given conditions is:")
    sympy.pretty_print(M_N_simplified)
    print("\nIts eigenvalues are:", list(eigenvalues.keys()))
    print("\nThe eigenvalue not proportional to M1, M2, or mu is:")
    print(target_eigenvalue)

solve_neutralino_mass()