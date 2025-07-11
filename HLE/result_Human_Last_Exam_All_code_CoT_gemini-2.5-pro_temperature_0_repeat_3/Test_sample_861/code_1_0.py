import sympy

def solve_difference_equation():
    """
    Solves the given difference equation and calculates the final expression.
    """
    # Define symbols for the variables in the equation
    n = sympy.Symbol('n')
    A, C = sympy.symbols('A C')

    # Step 1: Find the homogeneous solution by solving the characteristic equation
    # 8*r**2 - 6*r + 1 = 0
    r = sympy.Symbol('r')
    char_eq = 8*r**2 - 6*r + 1
    roots = sympy.solve(char_eq, r)

    # Step 2: Assign roots to B and D.
    # By convention, we assign the root with the larger magnitude to B.
    if abs(roots[0]) > abs(roots[1]):
        B = roots[0]
        D = roots[1]
    else:
        B = roots[1]
        D = roots[0]

    # Step 3: Find the particular solution E.
    # Assume y_p[n] = E (a constant). Substitute into the original equation:
    # 8*E - 6*E + E = 1
    E_sym = sympy.Symbol('E_sym')
    particular_eq = 8*E_sym - 6*E_sym + E_sym - 1
    E = sympy.solve(particular_eq, E_sym)[0]

    # Step 4: Use initial conditions to find A and C.
    # The general solution is y[n] = A*B**n + C*D**n + E
    y_n = A*B**n + C*D**n + E

    # Given initial conditions
    y0 = 1
    y_minus_1 = 2

    # Set up a system of two linear equations
    # For n=0: y[0] = A + C + E = 1
    # For n=-1: y[-1] = A/B + C/D + E = 2
    eq1 = sympy.Eq(y_n.subs(n, 0), y0)
    eq2 = sympy.Eq(y_n.subs(n, -1), y_minus_1)

    # Solve the system for A and C
    solution = sympy.solve((eq1, eq2), (A, C))
    A_val = solution[A]
    C_val = solution[C]

    # Print the components of the final equation as requested
    print("The closed form solution is y[n] = A * B^n + C * D^n + E")
    print("The determined coefficients are:")
    print(f"A = {A_val}")
    print(f"B = {B}")
    print(f"C = {C_val}")
    print(f"D = {D}")
    print(f"E = {E}")
    print("\nSo the final equation is:")
    print(f"y[n] = ({A_val}) * ({B})^n + ({C_val}) * ({D})^n + {E}")

    # Step 5: Calculate the value of the final expression
    result = E/A_val + (D*C_val)/B
    
    print("\nCalculating the expression E/A + (D*C)/B:")
    print(f"E/A + (D*C)/B = ({E})/({A_val}) + (({D})*({C_val}))/({B}) = {result}")
    print(f"The numerical value is: {float(result)}")

if __name__ == '__main__':
    solve_difference_equation()