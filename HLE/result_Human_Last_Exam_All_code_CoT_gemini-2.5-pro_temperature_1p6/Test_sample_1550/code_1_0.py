import sympy

def count_real_intersections():
    """
    This function calculates the number of real intersection points
    for the two given quadric equations.
    """
    # Define x and y as symbolic variables
    x, y = sympy.symbols('x y')

    # Define the first quadric equation
    # 164*x^2 - 216*x*y + 72*y^2 - 16*x + 31 = 0
    eq1 = sympy.Eq(164*x**2 - 216*x*y + 72*y**2 - 16*x + 31, 0)

    # Define the second quadric equation
    # 864*x^2 - 1056*x*y + 324*y^2 - 560*x + 324*y + 149 = 0
    eq2 = sympy.Eq(864*x**2 - 1056*x*y + 324*y**2 - 560*x + 324*y + 149, 0)

    # Solve the system of equations for x and y
    solutions = sympy.solve([eq1, eq2], [x, y])

    # Count the number of solutions where both x and y are real numbers
    real_solution_count = 0
    for sol in solutions:
        # A solution is a tuple (x_val, y_val).
        # We check if both values are real using the .is_real property.
        if sol[0].is_real and sol[1].is_real:
            real_solution_count += 1
    
    # Per the instructions, we output the single number representing the final count.
    # The 'final equation' is conceptually 'count = N', so we output 'N'.
    print(real_solution_count)

count_real_intersections()