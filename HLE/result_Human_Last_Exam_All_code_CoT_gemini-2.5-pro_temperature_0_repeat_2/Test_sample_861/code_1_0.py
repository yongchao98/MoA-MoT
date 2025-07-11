import sympy

def solve_difference_equation():
    """
    Solves the given linear constant-coefficient difference equation
    and calculates the required expression.
    """
    # Define symbols for the unknown coefficients
    A, C = sympy.symbols('A C')

    # --- Step 1: Find the particular solution E ---
    # The equation is 8y[n] - 6y[n-1] + y[n-2] = 1.
    # Assume a particular solution yp[n] = E (a constant).
    # Substitute into the equation: 8*E - 6*E + E = 1 => 3*E = 1
    E = sympy.Rational(1, 3)

    # --- Step 2: Find the homogeneous solution roots B and D ---
    # The characteristic equation is 8r^2 - 6r + 1 = 0.
    r = sympy.symbols('r')
    char_eq = 8*r**2 - 6*r + 1
    roots = sympy.solve(char_eq, r)

    # By convention, assign the larger root to B and the smaller to D.
    B = max(roots)
    D = min(roots)

    # --- Step 3 & 4: Use initial conditions to find A and C ---
    # The total solution is y[n] = A*B**n + C*D**n + E.
    # Initial conditions: y[0] = 1 and y[-1] = 2.
    # For n=0: A*B**0 + C*D**0 + E = 1  => A + C + E = 1
    # For n=-1: A*B**-1 + C*D**-1 + E = 2
    eq1 = sympy.Eq(A + C + E, 1)
    eq2 = sympy.Eq(A * B**-1 + C * D**-1 + E, 2)

    # Solve the system of equations for A and C
    solution = sympy.solve((eq1, eq2), (A, C))
    A_val = solution[A]
    C_val = solution[C]

    # --- Step 5: Calculate the final expression ---
    # The expression to calculate is E/A + (D*C)/B
    final_expression_value = E / A_val + (D * C_val) / B

    # --- Print the results ---
    print("Solving the difference equation: 8y[n] - 6y[n-1] + y[n-2] = 1")
    print(f"With initial conditions: y[0] = 1, y[-1] = 2")
    print("-" * 40)
    print("The constants for the solution y[n] = A*(B)^n + C*(D)^n + E are:")
    print(f"A = {A_val}")
    print(f"B = {B} (larger root of characteristic equation)")
    print(f"C = {C_val}")
    print(f"D = {D} (smaller root of characteristic equation)")
    print(f"E = {E} (from the particular solution)")
    print("-" * 40)
    print("The final closed-form solution is:")
    print(f"y[n] = ({A_val}) * ({B})^n + ({C_val}) * ({D})^n + ({E})")
    print("-" * 40)
    print("Calculating the expression E/A + (D*C)/B:")
    print(f"= ({E})/({A_val}) + (({D})*({C_val}))/({B})")
    print(f"= {E/A_val} + {D*C_val/B}")
    print(f"= {final_expression_value}")
    print("-" * 40)
    print(f"The final numerical value is: {float(final_expression_value)}")

if __name__ == '__main__':
    solve_difference_equation()