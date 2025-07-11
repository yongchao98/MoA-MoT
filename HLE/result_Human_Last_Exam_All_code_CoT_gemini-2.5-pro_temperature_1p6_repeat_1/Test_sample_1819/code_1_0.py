import sympy

def solve_flux():
    """
    This function calculates the energy flow through the yellow sides of the pyramid
    using the Divergence Theorem.
    """
    # Define the variables for integration
    x, y, z = sympy.symbols('x y z')

    # The vector field is F = (3x^3*y^2*z, 3x^2*y^3, z).
    # The divergence of F is div_F = 9*x^2*y^2*z + 9*x^2*y^2 + 1.
    div_F = 9*x**2*y**2*(z + 1) + 1

    # The pyramid's volume is defined by the integration bounds.
    # At a given height z, the cross-section is a square where x and y vary
    # between -s and s, with s = 1 - z/4.
    s = 1 - z/4

    # Integrate the divergence over the volume of the pyramid.
    # We integrate with respect to x first, then y, then z.
    
    # Integrate with respect to x:
    integral_over_x = sympy.integrate(div_F, (x, -s, s))
    
    # Integrate with respect to y:
    integral_over_xy = sympy.integrate(integral_over_x, (y, -s, s))
    
    # Integrate with respect to z to get the total flux through the closed surface:
    total_flux = sympy.integrate(integral_over_xy, (z, 0, 4))
    
    # The flux through the base (z=0) is 0.
    # So, total_flux is the flux through the four side faces.
    # The problem asks for the flux through the two yellow sides.
    # By symmetry of the problem's statement, this must be half of the total side flux.
    flux_yellow_sides = total_flux / 2

    # Output the steps of the final calculation
    print("The divergence of the field F is 9*x^2*y^2*(z + 1) + 1.")
    print("Integrating the divergence over the pyramid's volume gives the total flux.")
    print(f"The total flux through all four sides of the pyramid is: {total_flux}")
    print("The flux through the two yellow sides is half of this total.")
    print("Final Calculation:")
    print(f"{total_flux} / 2 = {flux_yellow_sides}")

# Execute the function to find the answer
solve_flux()