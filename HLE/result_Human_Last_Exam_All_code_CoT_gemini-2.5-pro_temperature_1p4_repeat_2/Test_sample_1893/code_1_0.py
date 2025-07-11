import sympy

def solve_neutralino_eigenvalue():
    """
    Solves for the special eigenvalue of the neutralino mass matrix based on the problem description.
    """
    # Define symbols for the parameters
    M1, M2, mu, beta, theta_W, MZ = sympy.symbols('M1 M2 mu beta theta_W MZ', real=True)
    sw = sympy.sin(theta_W)
    cw = sympy.cos(theta_W)

    # Define the neutralino mass matrix from the problem description
    M_11 = M1 * cw**2 + M2 * sw**2
    M_12 = (M2 - M1) * sw * cw
    M_22 = M1 * sw**2 + M2 * cw**2
    M_33 = mu * sympy.sin(2*beta)
    M_34 = -mu * sympy.cos(2*beta)
    M_44 = -mu * sympy.sin(2*beta)

    M_neutralino = sympy.Matrix([
        [M_11, M_12, 0, 0],
        [M_12, M_22, MZ, 0],
        [0, MZ, M_33, M_34],
        [0, 0, M_34, M_44]
    ])

    # Step 1: Apply the condition of "pure states"
    # This implies M1 = M2 and cos(2*beta) = 0.
    # For 0 <= beta <= pi/2, beta must be pi/4.
    M_simplified = M_neutralino.subs([(M2, M1), (beta, sympy.pi/4)])

    # Step 2: Look for an eigenvalue not proportional to M1, M2, or mu.
    # This is found in the special case where M1 = mu = 0.
    M_final = M_simplified.subs([(M1, 0), (mu, 0)])

    # Step 3: Calculate the eigenvalues of the final matrix
    eigenvals = M_final.eigenvals() # Returns a dictionary {eigenvalue: multiplicity}

    # The characteristic equation for the non-zero part is lambda^2 - MZ^2 = 0
    # Its solutions are lambda = +MZ and lambda = -MZ.
    
    # We are asked to compute the eigenvalue. We will use the known value for MZ.
    MZ_value = 91.1876  # Mass of Z boson in GeV/c^2

    print("The condition that tilde(gamma) and tilde(H)_b are pure states requires M1=M2 and beta=pi/4.")
    print("To find an eigenvalue not proportional to M1, M2, or mu, we consider the case where M1 = mu = 0.")
    print("In this case, the mass matrix simplifies significantly.")
    print("\nThe final matrix form is:")
    sympy.pprint(M_final)
    
    print(f"\nThe eigenvalues of this matrix are {list(eigenvals.keys())}.")
    
    # The final equation for the non-zero eigenvalues is lambda^2 = MZ^2
    # So lambda = +/- MZ
    
    non_zero_eigenvalues = [val for val in eigenvals.keys() if val != 0]

    # As requested, outputting the numbers in the final equation lambda = +/- MZ_value
    print("\nThe final equation for the non-zero eigenvalues (lambda) gives:")
    for val in non_zero_eigenvalues:
        # substitute MZ with its value for printing the final number
        num_val = val.subs(MZ, MZ_value)
        print(f"lambda = {num_val}")
    
    # The problem asks for "the" eigenvalue, which is typically the positive mass value.
    positive_eigenvalue = MZ_value
    print(f"\nThe computed positive eigenvalue is {positive_eigenvalue}.")

solve_neutralino_eigenvalue()