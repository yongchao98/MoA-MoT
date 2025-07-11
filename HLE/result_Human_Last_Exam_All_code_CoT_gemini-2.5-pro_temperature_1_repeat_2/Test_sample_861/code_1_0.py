import sympy

def solve_difference_equation():
    """
    Solves the given difference equation and calculates the final expression.
    """
    # Step 1: Find the homogeneous solution by solving the characteristic equation
    # The characteristic equation is 8r^2 - 6r + 1 = 0.
    r = sympy.Symbol('r')
    char_eq = 8*r**2 - 6*r + 1
    roots = sympy.solve(char_eq, r)
    
    # Assign the roots to B and D, with B being the larger root by convention.
    B, D = sorted(roots, reverse=True)

    # Step 2: Find the particular solution
    # Assume y_p[n] = E (a constant). Substitute into the original equation:
    # 8*E - 6*E + E = 1  =>  3*E = 1
    E = sympy.Rational(1, 3)

    # Step 3 & 4: Use initial conditions to find A and C
    # The general solution is y[n] = A*B**n + C*D**n + E.
    # We set up a system of linear equations for A and C.
    A = sympy.Symbol('A')
    C = sympy.Symbol('C')

    # From y[0] = 1:
    # A*(B**0) + C*(D**0) + E = 1  =>  A + C + E = 1
    eq1 = sympy.Eq(A + C + E, 1)

    # From y[-1] = 2:
    # A*(B**-1) + C*(D**-1) + E = 2  =>  A/B + C/D + E = 2
    eq2 = sympy.Eq(A/B + C/D + E, 2)

    # Step 5: Solve the system for A and C
    solution = sympy.solve((eq1, eq2), (A, C))
    A_val = solution[A]
    C_val = solution[C]

    print("--- Solving the Difference Equation ---")
    print("The equation is: 8y[n] - 6y[n-1] + y[n-2] = 1")
    print("Initial conditions: y[0] = 1, y[-1] = 2")
    print("The solution is of the form: y[n] = A*(B)^n + C*(D)^n + E\n")

    print(f"1. The roots of the characteristic equation (8r^2 - 6r + 1 = 0) are {B} and {D}.")
    print(f"   So, B = {B} and D = {D}.")
    print(f"2. The particular solution is a constant E = {E}.")
    print(f"3. Solving for A and C using initial conditions yields A = {A_val} and C = {C_val}.\n")

    # Print the final equation as requested
    print("--- Final Closed-Form Solution ---")
    print(f"y[n] = ({A_val}) * ({B})^n + ({C_val}) * ({D})^n + ({E})\n")

    # Step 6: Calculate the final expression
    final_expression_value = E/A_val + (D*C_val)/B
    
    print("--- Calculation of the Final Expression ---")
    print(f"The expression to calculate is E/A + (D*C)/B")
    print(f"= ({E}) / ({A_val}) + (({D}) * ({C_val})) / ({B})")
    print(f"= {E/A_val} + {D*C_val/B}")
    print(f"= {final_expression_value}\n")
    
    # Final answer in float format
    print(f"The final numerical answer is: {float(final_expression_value)}")


if __name__ == '__main__':
    solve_difference_equation()
<<<0.75>>>