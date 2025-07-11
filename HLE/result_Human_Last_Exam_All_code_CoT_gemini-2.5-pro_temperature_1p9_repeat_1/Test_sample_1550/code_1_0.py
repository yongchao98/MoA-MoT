import sympy

def solve_quadric_intersections():
    """
    This function finds the number of real intersection points for two given quadric equations.
    """
    # Define the variables
    x, y = sympy.symbols('x y')

    # Define the coefficients and constants for the first equation
    c1_x2, c1_xy, c1_y2, c1_x, c1_c = 164, -216, 72, -16, 31
    # Define the first quadric equation
    eq1 = c1_x2*x**2 + c1_xy*x*y + c1_y2*y**2 + c1_x*x + c1_c
    
    # Define the coefficients and constants for the second equation
    c2_x2, c2_xy, c2_y2, c2_x, c2_y, c2_c = 864, -1056, 324, -560, 324, 149
    # Define the second quadric equation
    eq2 = c2_x2*x**2 + c2_xy*x*y + c2_y2*y**2 + c2_x*x + c2_y*y + c2_c

    # Print the equations with their numbers
    print(f"The first quadric is given by the equation: {c1_x2}*x^2 + ({c1_xy})*x*y + {c1_y2}*y^2 + ({c1_x})*x + {c1_c} = 0")
    print(f"The second quadric is given by the equation: {c2_x2}*x^2 + ({c2_xy})*x*y + {c2_y2}*y^2 + ({c2_x})*x + {c2_y}*y + {c2_c} = 0")

    # Solve the system of equations
    try:
        solutions = sympy.solve([eq1, eq2], [x, y], dict=False)
    except Exception as e:
        print(f"An error occurred during solving: {e}")
        return

    # Count the number of real solutions
    real_solutions_count = 0
    if solutions:
        for sol in solutions:
            # Check if both components of the solution are real numbers.
            # is_real is a property of sympy expressions
            if sol[0].is_real and sol[1].is_real:
                real_solutions_count += 1
    
    print(f"\nThe number of real intersection points is: {real_solutions_count}")

solve_quadric_intersections()
<<<4>>>