import sympy
from sympy import integrate, symbols, Rational, pretty

def solve_flux():
    """
    Calculates the flux of a vector field F through the yellow sides of a pyramid.
    """
    # Define the symbols
    x, y, z = symbols('x y z')

    # The vector field is F = (3x^3*y^2*z, 3x^2*y^3, z)
    # We are calculating the flux through the two yellow sides.
    # Let's consider the side S1 in the x>0 region, defined by the plane 4x + z = 4.
    # We parameterize this surface by projecting it onto the yz-plane.
    # x = 1 - z/4
    # The integration domain in the yz-plane is a triangle with z from 0 to 4,
    # and for each z, y goes from -(1-z/4) to (1-z/4).
    # The outward normal vector element dS is (1, 0, 1/4) dy dz.

    # The dot product F . dS is (3*x**3*y**2*z) * 1 + z * (1/4)
    # Substitute x = 1 - z/4
    integrand_part_1 = 3 * (1 - z/4)**3 * y**2 * z
    integrand_part_2 = z / 4
    
    integrand = integrand_part_1 + integrand_part_2

    # Define integration limits for y
    y_lower = -(1 - z/4)
    y_upper = 1 - z/4

    # Integrate with respect to y first, then z
    # We can integrate the two parts separately for clarity
    
    # Integral of the first part
    integral_1 = integrate(integrand_part_1, (y, y_lower, y_upper), (z, 0, 4))
    
    # Integral of the second part
    integral_2 = integrate(integrand_part_2, (y, y_lower, y_upper), (z, 0, 4))

    # The flux through one yellow side (S1) is the sum of these two integrals.
    flux_one_side = integral_1 + integral_2

    # The other yellow side (S2, at x<0) is symmetric, and due to the form of F and dS,
    # the flux through it is identical to the flux through S1.
    # F.dS_2 = -3*x**3*y**2*z + z/4. With x=-(1-z/4), this becomes
    # -3*(-(1-z/4))**3*y**2*z + z/4 = 3*(1-z/4)**3*y**2*z + z/4, which is the same integrand.
    # So, the total flux is 2 * flux_one_side.
    
    total_flux = 2 * flux_one_side

    # Output the results in a clear equation format
    flux1_str = f"{integral_1.p} / {integral_1.q}" if isinstance(integral_1, Rational) else str(integral_1)
    flux2_str = f"{integral_2.p} / {integral_2.q}" if isinstance(integral_2, Rational) else str(integral_2)
    flux_one_side_str = f"{flux_one_side.p} / {flux_one_side.q}" if isinstance(flux_one_side, Rational) else str(flux_one_side)
    total_flux_str = f"{total_flux.p} / {total_flux.q}" if isinstance(total_flux, Rational) else str(total_flux)

    print("The flux through one yellow side is the sum of two integrals:")
    print(f"Integral 1 = {flux1_str}")
    print(f"Integral 2 = {flux2_str}")
    print(f"Flux through one side = {flux1_str} + {flux2_str} = {flux_one_side_str}\n")
    
    print("The flux through the other yellow side is identical.")
    print("The total energy flow through both yellow sides is the sum of the flux through each side.")
    print("Final Equation:")
    print(f"{flux_one_side_str} + {flux_one_side_str} = {total_flux_str}")
    
    # The final answer in the required format
    print(f"\n<<<{total_flux_str}>>>")

solve_flux()