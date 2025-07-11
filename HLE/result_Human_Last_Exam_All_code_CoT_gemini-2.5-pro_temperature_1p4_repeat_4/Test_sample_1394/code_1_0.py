import sympy

def display_solution():
    """
    This function defines and prints the general solution for the given ODE.
    The final equation is y^2 = C*x^2 + C^2 + C + 9, where C is an arbitrary constant.
    The numbers in the equation are explicitly part of the symbolic representation.
    """
    # Define symbols for the variables and the arbitrary constant
    x, y, C = sympy.symbols('x y C')

    # Construct the equation for the general solution
    # The left-hand side is y squared
    lhs = y**2
    # The right-hand side is a polynomial in x and C
    rhs = C * x**2 + C**2 + C + 9

    solution_equation = sympy.Eq(lhs, rhs)

    # Print the resulting symbolic equation
    print("The general solution is:")
    print(solution_equation)

if __name__ == '__main__':
    display_solution()