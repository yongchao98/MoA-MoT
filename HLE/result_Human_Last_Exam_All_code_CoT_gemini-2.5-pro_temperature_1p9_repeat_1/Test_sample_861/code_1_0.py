import sympy

def solve_difference_equation_and_calculate():
    """
    This script solves the given linear difference equation and calculates the specified expression.
    The steps are:
    1. Find the particular solution (constant E).
    2. Find the homogeneous solution by solving the characteristic equation for its roots (B and D).
    3. Use the initial conditions to set up and solve a system of linear equations for coefficients A and C.
    4. Substitute the found values into the expression E/A + (D*C)/B and calculate the result.
    """

    # --- Step 1: Find the particular solution E ---
    # The equation is 8y[n] - 6y[n-1] + y[n-2] = 1.
    # We assume a constant particular solution y_p[n] = E.
    # Substituting into the equation gives: 8*E - 6*E + E = 1 => 3*E = 1.
    E = sympy.Rational(1, 3)

    # --- Step 2: Find the homogeneous solution roots B and D ---
    # The homogeneous equation is 8y_h[n] - 6y_h[n-1] + y_h[n-2] = 0.
    # This leads to the characteristic equation: 8*r**2 - 6*r + 1 = 0.
    r = sympy.Symbol('r')
    char_eq = 8*r**2 - 6*r + 1
    roots = sympy.solve(char_eq, r)

    # The solution form is y[n] = A*(B)**n + C*(D)**n + E.
    # By convention, we assign the larger root to B and the smaller to D.
    B = max(roots)
    D = min(roots)

    # --- Step 3: Find the coefficients A and C using initial conditions ---
    # Initial conditions are given as y[0] = 1 and y[-1] = 2.
    # The general solution is y[n] = A*(B**n) + C*(D**n) + E.
    A_sym, C_sym = sympy.symbols('A C')

    # Using y[0] = 1: A*(B**0) + C*(D**0) + E = 1 => A + C + E = 1
    eq1 = sympy.Eq(A_sym + C_sym + E, 1)

    # Using y[-1] = 2: A*(B**-1) + C*(D**-1) + E = 2 => A/B + C/D + E = 2
    eq2 = sympy.Eq(A_sym/B + C_sym/D + E, 2)

    # Solve the system of two linear equations for A and C
    solution = sympy.solve((eq1, eq2), (A_sym, C_sym))
    A = solution[A_sym]
    C = solution[C_sym]

    # --- Step 4: Output the constants and calculate the final expression ---
    print("The solution to the difference equation is y[n] = A*(B)^n + C*(D)^n + E")
    print("The calculated values for the constants are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}")
    print("-" * 30)

    # As requested, output the final equation with the numbers
    print("The final closed-form solution is:")
    print(f"y[n] = ({A}) * ({B})^n + ({C}) * ({D})^n + ({E})")
    print("-" * 30)
    
    # Calculate the final expression E/A + (D*C)/B
    expression_val = E/A + (D*C)/B

    print("Calculating the value of the expression E/A + (D*C)/B:")
    term1 = E/A
    term2 = (D*C)/B
    print(f"The first term E/A is ({E})/({A}) = {term1}")
    print(f"The second term (D*C)/B is (({D})*({C}))/({B}) = {term2}")
    print(f"The final result is {term1} + {term2} = {expression_val}")

if __name__ == "__main__":
    solve_difference_equation_and_calculate()