import sympy

def solve_neutralino_mass():
    """
    This function solves for the neutralino mass eigenvalue under specific conditions.
    
    The steps are:
    1. Define the neutralino mass matrix symbolically.
    2. Apply physical conditions for "dynamic enhancement" which cause the matrix to become block-diagonal.
       - Condition 1: M1 = M2. Derived from setting the mixing term (M2 - M1)*sin(theta_W)*cos(theta_W) to zero.
       - Condition 2: beta = pi/4. Derived from setting the mixing term -mu*cos(2*beta) to zero.
    3. Apply the condition that the eigenvalue is not proportional to M1 or mu, which implies a scenario where M1=mu=0.
    4. Formulate and solve the characteristic equation for the resulting simplified matrix.
    5. The non-zero solutions will be the eigenvalues that depend only on M_Z.
    """
    
    # Define symbols for the parameters
    M1, M2, mu, MZ, b, l = sympy.symbols('M1 M2 mu M_Z beta lambda')

    # Define the original matrix from the problem description
    sWcW = sympy.Symbol('swcw') # A placeholder for sin(theta_W)*cos(theta_W)
    cW2 = sympy.Symbol('cw2')   # A placeholder for cos^2(theta_W)
    sW2 = sympy.Symbol('sw2')   # A placeholder for sin^2(theta_W)

    # Full matrix as given
    M_full = sympy.Matrix([
      [ M1*cW2 + M2*sW2, (M2 - M1)*sWcW,   0,                   0                   ],
      [ (M2 - M1)*sWcW,  M1*sW2 + M2*cW2,  MZ,                  0                   ],
      [ 0,                 MZ,                 mu*sympy.sin(2*b), -mu*sympy.cos(2*b)    ],
      [ 0,                 0,                  -mu*sympy.cos(2*b),  -mu*sympy.sin(2*b)    ]
    ])

    # 1. Apply condition M1 = M2
    M_step1 = M_full.subs(M2, M1)
    
    # Since cW2+sW2 = 1, the diagonal terms M1*cW2 + M1*sW2 simplify to M1
    M_step1[0,0] = M1
    M_step1[1,1] = M1

    # 2. Apply condition beta = pi/4, so cos(2*beta) = 0 and sin(2*beta) = 1
    M_step2 = M_step1.subs(b, sympy.pi/4)

    # 3. Apply condition for independence from M1, mu => M1=0, mu=0
    M_final = M_step2.subs({M1: 0, mu: 0})
    
    # 4. Calculate the characteristic equation: det(M_final - lambda*I) = 0
    char_poly = M_final.charpoly(l)
    
    # Get coefficients to display the full equation form
    coeffs = char_poly.all_coeffs()
    c4, c3, c2, c1, c0 = coeffs

    # The mass of the Z boson in GeV
    MZ_val = 91.1876
    
    c2_val = c2.subs(MZ, MZ_val)

    print("The final characteristic equation for the eigenvalue (λ) is:")
    print("a*λ^4 + b*λ^3 + c*λ^2 + d*λ + e = 0")
    print("\nWhere the coefficients are:")
    print(f"a = {c4}")
    print(f"b = {c3}")
    print(f"c = {c2_val}")
    print(f"d = {c1}")
    print(f"e = {c0}")
    
    # 5. Find the eigenvalues by solving the equation
    # The non-zero eigenvalues are the ones we seek
    eigenvals = sympy.solve(sympy.Eq(char_poly.as_expr(), 0), l)
    
    target_eigenvals = [val for val in eigenvals if val != 0]

    # The question asks for 'the' eigenvalue. Conventionally, particle masses are positive.
    positive_eigenvalue = max(target_eigenvals)
    result = positive_eigenvalue.subs(MZ, MZ_val)

    print(f"\nThe resulting eigenvalue that is not proportional to M1, M2, or mu is {result} GeV.")

solve_neutralino_mass()
<<<91.1876>>>