import sympy

def solve_flux():
    """
    Calculates the energy flow through the yellow sides of the pyramid.
    """
    # Define the symbolic variables
    x, z = sympy.symbols('x z')

    print("Step 1: Set up the integral for the total flux.")
    print("The total energy flow (flux) through the two yellow sides is given by the integral:")
    # Using sympy.pretty_print for a nicer formula layout
    equation_text = "Flux = (1/2) * Integral(Integral((3/16)*x**2*(4-z)**3 + z, x), z)"
    print(equation_text)
    print("The integration bounds are z from 0 to 4, and x from -(1 - z/4) to (1 - z/4).\n")

    # Define the integrand
    integrand = (sympy.S(3)/16) * x**2 * (4 - z)**3 + z

    print("Step 2: Calculate the inner integral with respect to x.")
    # Integrate the expression with respect to x
    inner_integral_result = sympy.integrate(integrand, (x, -(1 - z/4), (1 - z/4)))
    inner_integral_simplified = sympy.simplify(inner_integral_result)
    
    # We will show the two terms of the result separately for clarity
    term1_after_x_integration = sympy.simplify((sympy.S(3)/16) * (4-z)**3 * sympy.integrate(x**2, (x, -(1 - z/4), (1 - z/4))))
    term2_after_x_integration = sympy.simplify(z * sympy.integrate(1, (x, -(1 - z/4), (1 - z/4))))

    print(f"The integral of the first term ( (3/16)*x**2*(4-z)**3 ) with respect to x is: {term1_after_x_integration}")
    print(f"The integral of the second term ( z ) with respect to x is: {term2_after_x_integration}")
    print(f"The result of the inner integral is the sum: {inner_integral_simplified}\n")


    print("Step 3: Calculate the outer integral with respect to z.")
    # Integrate the result from the inner integral with respect to z
    # We integrate the two parts separately to show the numbers
    part1_integral_z = sympy.integrate(term1_after_x_integration, (z, 0, 4))
    part2_integral_z = sympy.integrate(term2_after_x_integration, (z, 0, 4))
    total_integral_value = part1_integral_z + part2_integral_z
    
    print(f"Integrating the first part ({term1_after_x_integration}) from z=0 to 4 gives: {part1_integral_z}")
    print(f"Integrating the second part ({term2_after_x_integration}) from z=0 to 4 gives: {part2_integral_z}")
    print(f"The total value of the double integral is {part1_integral_z} + {part2_integral_z} = {total_integral_value}\n")

    print("Step 4: Final Calculation for the flux.")
    # The total flux is 1/2 of the integral's value
    total_flux = sympy.S(1)/2 * total_integral_value
    print(f"The total flux is (1/2) times the value of the integral.")
    print(f"Total Flux = (1/2) * {total_integral_value} = {total_flux}")

solve_flux()