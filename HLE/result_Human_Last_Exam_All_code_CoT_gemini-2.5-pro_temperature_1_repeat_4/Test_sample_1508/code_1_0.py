import sympy

def demonstrate_linear_independence():
    """
    Constructs an example for n=3, s=2 and demonstrates the linear independence
    of the corresponding polynomials P_i(x).
    """
    # Setup parameters
    n = 3
    s = 2  # s > floor(n/2) = 1
    L = {0, 1}

    print(f"Let n = {n}, s = {s}. The condition s > floor(n/2) holds since {s} > {n//2}.")
    print(f"Let L = {L}.")

    # Define an ordered L-intersecting family F
    # F_1 = {1, 3}, F_2 = {2, 3}.
    # n=3 is in both, so r=2. |F_1|=|F_2|=2, so size ordering holds.
    # |F_1 intersect F_2| = |{3}| = 1, which is in L.
    F1 = {1, 3}
    F2 = {2, 3}
    v1 = [1, 0, 1]
    v2 = [0, 1, 1]
    
    print(f"Consider the ordered L-intersecting family F = {{ {F1}, {F2} }}.")
    
    # Define sympy variables
    x = sympy.symbols('x_1:{}'.format(n + 1))
    
    # Define the polynomials P_i(x)
    # For F1, |F1|=2. l_k < 2 are l=0 and l=1.
    # P1(x) = (<x,v1> - 0) * (<x,v1> - 1)
    v1_dot_x = sum(vi * xi for vi, xi in zip(v1, x))
    P1 = v1_dot_x * (v1_dot_x - 1)
    
    # For F2, |F2|=2. l_k < 2 are l=0 and l=1.
    # P2(x) = (<x,v2> - 0) * (<x,v2> - 1)
    v2_dot_x = sum(vi * xi for vi, xi in zip(v2, x))
    P2 = v2_dot_x * (v2_dot_x - 1)
    
    print("\nThe corresponding polynomials are:")
    print(f"P_1(x) = {sympy.expand(P1)}")
    print(f"P_2(x) = {sympy.expand(P2)}")
    
    # Check for linear dependence
    # c1*P1 + c2*P2 = 0
    c1, c2 = sympy.symbols('c1 c2')
    equation = sympy.expand(c1 * P1 + c2 * P2)
    
    print("\nTo check for linear dependence, we solve the equation c1*P1 + c2*P2 = 0:")
    print(f"{equation} = 0")
    
    # This polynomial must be zero for all x, so all its coefficients must be zero.
    # We can collect coefficients for each monomial.
    coeffs = sympy.Poly(equation, x).coeffs()
    
    # Or, more systematically, solve for c1, c2
    solution = sympy.solve(equation, (c1, c2), dict=True)
    
    print("\nSetting the coefficients of the monomials to zero gives a system of linear equations in c1, c2.")
    print("For example:")
    print(f"Coefficient of x_1**2: {equation.coeff(x[0]**2)} = 0  => c1 = 0")
    print(f"Coefficient of x_2**2: {equation.coeff(x[1]**2)} = 0  => c2 = 0")

    if not solution or all(s[c] == 0 for s in solution for c in (c1, c2)):
        print("\nThe only solution is c1 = 0, c2 = 0.")
        print("This shows that the polynomials are linearly independent.")
    else:
        print("\nFound a non-trivial solution, the polynomials are linearly dependent.")
        print(solution)

demonstrate_linear_independence()