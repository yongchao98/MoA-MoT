import sympy

def solve_quadric_system():
    """
    This function finds the number of real intersection points for the given system of two quadric equations.
    """
    # Define the variables
    x, y = sympy.symbols('x y', real=True)

    # Define the coefficients of the first equation: 164*x^2 - 216*x*y + 72*y^2 - 16*x + 31 = 0
    a1, b1, c1, d1, e1, f1 = 164, -216, 72, -16, 0, 31
    # Define the coefficients of the second equation: 864*x^2 - 1056*x*y + 324*y^2 - 560*x + 324*y + 149 = 0
    a2, b2, c2, d2, e2, f2 = 864, -1056, 324, -560, 324, 149

    # Create the equations using the specified numbers
    eq1 = a1*x**2 + b1*x*y + c1*y**2 + d1*x + e1*y + f1
    eq2 = a2*x**2 + b2*x*y + c2*y**2 + d2*x + e2*y + f2

    print("Solving the system of equations:")
    print(f"Equation 1: {a1}*x**2 + ({b1})*x*y + {c1}*y**2 + ({d1})*x + {e1}*y + {f1} = 0")
    print(f"Equation 2: {a2}*x**2 + ({b2})*x*y + {c2}*y**2 + ({d2})*x + {e2}*y + {f2} = 0")
    print("-" * 20)
    
    # Solve the system of equations
    # Using 'dict=True' returns a list of dictionaries, which is easier to parse.
    solutions = sympy.solve([eq1, eq2], [x, y], dict=True)

    # Filter for solutions where both x and y are real numbers and count them.
    # The 'real=True' flag in sympy.symbols helps, but it's good practice to verify.
    real_solutions_count = 0
    for sol in solutions:
        # Check if the solution is fully real. A solution might be empty or contain non-real components.
        if sol and all(val.is_real for val in sol.values()):
            real_solutions_count += 1
            
    print(f"Number of real intersection points: {real_solutions_count}")

if __name__ == '__main__':
    solve_quadric_system()