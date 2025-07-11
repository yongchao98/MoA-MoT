import sympy

def solve_neutralino_mass():
    """
    Solves for the specific eigenvalue of the neutralino mass matrix based on the problem's conditions.
    """
    # Define symbols for the parameters in the mass matrix
    M1, M2, mu, beta, MZ = sympy.symbols('M1 M2 mu beta MZ', real=True)
    sw, cw = sympy.symbols('sw cw', real=True, positive=True)
    lambda_sym = sympy.symbols('lambda')

    # Define the full Neutralino mass matrix symbolically
    M_N = sympy.Matrix([
        [M1*cw**2 + M2*sw**2, (M2 - M1)*sw*cw, 0, 0],
        [(M2 - M1)*sw*cw, M1*sw**2 + M2*cw**2, MZ, 0],
        [0, MZ, mu*sympy.sin(2*beta), -mu*sympy.cos(2*beta)],
        [0, 0, -mu*sympy.cos(2*beta), -mu*sympy.sin(2*beta)]
    ])

    print("Step 1: Apply conditions for dynamic enhancement (M1=M2, beta=pi/4).")
    # Substitute M2=M1 and beta=pi/4 (so cos(2*beta)=0, sin(2*beta)=1)
    M_N_simplified = M_N.subs({M2: M1, sympy.cos(2*beta): 0, sympy.sin(2*beta): 1})
    print("Simplified Mass Matrix:")
    sympy.pprint(M_N_simplified)
    print("-" * 40)

    print("Step 2: Find the eigenvalues. The matrix is block-diagonal.")
    print("The eigenvalues are M1, -mu, and the eigenvalues of the central 2x2 block.")
    central_block = M_N_simplified[1:3, 1:3]
    print("Central 2x2 Block:")
    sympy.pprint(central_block)
    print("-" * 40)

    print("Step 3: Apply the condition that the eigenvalue is independent of M1 and mu.")
    print("This requires setting M1=0 and mu=0 in the central block.")
    final_block = central_block.subs({M1: 0, mu: 0})
    print("Final form of the central block:")
    sympy.pprint(final_block)
    print("-" * 40)
    
    print("Step 4: Compute the characteristic equation for this final block.")
    # The characteristic equation is det(M - lambda*I) = 0
    char_eq_expr = (final_block - lambda_sym * sympy.eye(2)).det()
    final_equation = sympy.Eq(char_eq_expr, 0)
    
    print("The final equation for the non-zero eigenvalues is:")
    sympy.pprint(final_equation)
    print("\nIn the format a*lambda^2 + b*lambda + c = 0, the equation is:")
    print(f"(1) * lambda**2 + (0) * lambda + (-1)*MZ**2 = 0")
    print("The numbers in the final equation are 1, 0, and -1.")
    print("-" * 40)

    print("Step 5: Solve the equation to find the eigenvalues.")
    eigenvalues = sympy.solve(final_equation, lambda_sym)
    
    print(f"The solutions are lambda = {eigenvalues[0]} and lambda = {eigenvalues[1]}.")
    print("The other two eigenvalues of the full 4x4 matrix are 0.")
    print("\nThe question asks for 'the' eigenvalue not proportional to M1, M2, or mu.")
    print("We provide the positive value, which corresponds to a physical mass.")
    
    # The positive eigenvalue is MZ
    final_answer = eigenvalues[1]
    
    print("\nThe computed eigenvalue is:")
    print(final_answer)

solve_neutralino_mass()