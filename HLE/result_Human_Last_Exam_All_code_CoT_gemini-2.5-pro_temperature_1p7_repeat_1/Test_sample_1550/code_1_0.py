import sympy

def solve_quadric_system():
    """
    This function finds the number of real intersection points for the two given quadric equations.
    """
    # Define the variables x and y, specifying they must be real.
    x, y = sympy.symbols('x y', real=True)

    # First quadric equation: 164*x^2 - 216*x*y + 72*y^2 - 16*x + 31 = 0
    # Store the equation as an expression. Note that we don't need the "= 0" part for the solver.
    c1 = [164, -216, 72, -16, 0, 31]
    eq1 = c1[0]*x**2 + c1[1]*x*y + c1[2]*y**2 + c1[3]*x + c1[4]*y + c1[5]
    eq1_str = f"{c1[0]}*x^2 + ({c1[1]})*x*y + {c1[2]}*y^2 + ({c1[3]})*x + {c1[5]} = 0"

    # Second quadric equation: 864*x^2 - 1056*x*y + 324*y^2 - 560*x + 324*y + 149 = 0
    # Store the equation as an expression.
    c2 = [864, -1056, 324, -560, 324, 149]
    eq2 = c2[0]*x**2 + c2[1]*x*y + c2[2]*y**2 + c2[3]*x + c2[4]*y + c2[5]
    eq2_str = f"{c2[0]}*x^2 + ({c2[1]})*x*y + {c2[2]}*y^2 + ({c2[3]})*x + {c2[4]}*y + {c2[5]} = 0"

    # Use nonlinsolve to find the set of real solutions for the system of equations.
    # nonlinsolve is suitable for systems of polynomial equations and returns a set of solutions.
    solutions = sympy.nonlinsolve([eq1, eq2], [x, y])

    # The number of real solutions is the number of elements in the resulting set.
    num_real_solutions = len(solutions)

    # Print the problem statement and the final answer.
    print("Finding the number of real intersection points for the plane quadrics:")
    print(f"Equation 1: {eq1_str}")
    print(f"Equation 2: {eq2_str}")
    print(f"\nThe number of real intersection points is: {num_real_solutions}")

if __name__ == '__main__':
    solve_quadric_system()