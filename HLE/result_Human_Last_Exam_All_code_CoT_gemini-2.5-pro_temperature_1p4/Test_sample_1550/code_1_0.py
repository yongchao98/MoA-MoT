import sympy

def count_real_intersections():
    """
    Solves a system of two quadric equations and counts the number of real solutions.
    """
    # Define the variables as symbolic objects
    x, y = sympy.symbols('x y')

    # Define the coefficients of the first quadric equation:
    # 164*x^2 - 216*x*y + 72*y^2 - 16*x + 31 = 0
    c1 = {'A': 164, 'B': -216, 'C': 72, 'D': -16, 'E': 0, 'F': 31}
    eq1 = c1['A']*x**2 + c1['B']*x*y + c1['C']*y**2 + c1['D']*x + c1['E']*y + c1['F']

    # Define the coefficients of the second quadric equation:
    # 864*x^2 - 1056*x*y + 324*y^2 - 560*x + 324*y + 149 = 0
    c2 = {'A': 864, 'B': -1056, 'C': 324, 'D': -560, 'E': 324, 'F': 149}
    eq2 = c2['A']*x**2 + c2['B']*x*y + c2['C']*y**2 + c2['D']*x + c2['E']*y + c2['F']

    # Print the original equations with their coefficients
    print("Finding the number of real intersection points for the quadrics:")
    print(f"Equation 1: {c1['A']}*x^2 + ({c1['B']})*x*y + {c1['C']}*y^2 + ({c1['D']})*x + {c1['E']}*y + {c1['F']} = 0")
    print(f"Equation 2: {c2['A']}*x^2 + ({c2['B']})*x*y + {c2['C']}*y^2 + ({c2['D']})*x + {c2['E']}*y + {c2['F']} = 0")
    print("\n" + "="*50)

    # Solve the system of equations.
    # The `solve` function finds all complex solutions.
    solutions = sympy.solve([eq1, eq2], [x, y], dict=True)

    # Filter for real solutions and count them.
    # A solution is real if both its x and y components are real numbers.
    # The `is_real` attribute of a SymPy expression is used for this check.
    real_solutions_count = 0
    if solutions:
        for sol in solutions:
            # Check if both x and y are real numbers.
            # `sol[x]` is the x-value of the solution, `sol[y]` is the y-value.
            if sol[x].is_real and sol[y].is_real:
                real_solutions_count += 1
    
    # Print the final result. If the solutions list is empty, there are no intersections.
    if not solutions:
        print("The system of equations has no solutions in the complex domain.")
        print("Therefore, there are 0 real intersection points.")
    else:
        print(f"The system of equations has {len(solutions)} solution(s) in the complex domain.")
        print(f"Out of these, there are {real_solutions_count} real intersection point(s).")

if __name__ == '__main__':
    count_real_intersections()