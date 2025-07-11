import sympy

def calculate_energy_flow():
    """
    Calculates the energy flow (flux) of the vector field F through the two yellow sides of a square pyramid.
    """
    # Define symbolic variables
    y, z = sympy.symbols('y z')

    # The vector field is F = (3x^3y^2z, 3x^2y^3, z).
    # The pyramid has its base on the xy-plane (z=0), centered at the origin, with side length 2.
    # The apex is at (0,0,4). The four triangular faces connect the base edges to the apex.
    # Let's assume the front face (in the +x direction) and back face (in the -x direction) are yellow.

    # We need to calculate the flux integral integral(F . dS) over the two yellow surfaces.
    # Due to the symmetry of the pyramid and the properties of F, the flux through the front face
    # is identical to the flux through the back face. We will calculate the flux for one face and multiply by 2.

    # --- Calculation for the Front Face ---
    # The front face is a triangle with base vertices (1,-1,0) and (1,1,0), and apex (0,0,4).
    # The plane containing this face has the equation 4x + z = 4.
    # An outward normal vector is n = (4, 0, 1).
    # We parameterize the surface by projecting it onto the yz-plane.
    # The differential surface element is dS = n dy dz / |n . i| = (4, 0, 1) dy dz / 4.
    # On the surface, we substitute x = (4-z)/4 = 1 - z/4.
    
    # The integrand for the flux is (F . n) / 4.
    # F . n = (3x^3y^2z, 3x^2y^3, z) . (4, 0, 1) = 12x^3y^2z + z
    # After substituting x = 1 - z/4:
    # F . n = 12 * (1 - z/4)**3 * y**2 * z + z
    # The integrand is (F . n) / 4 = 3 * z * y**2 * (1 - z/4)**3 + z/4

    integrand = 3 * z * y**2 * (1 - z/4)**3 + z/4

    # The projection on the yz-plane is a triangle with vertices (y,z) = (-1,0), (1,0), (0,4).
    # The integration limits for z are [0, 4]. For a given z, y ranges from -(1-z/4) to (1-z/4).
    y_limit = 1 - z/4
    
    # Integrate with respect to y
    integral_vs_y = sympy.integrate(integrand, (y, -y_limit, y_limit))

    # The result of the first integration is the function we need to integrate over z.
    # This comes out to: 2*z*(1 - z/4)**6 + (z/2)*(1 - z/4)
    # This can be split into two parts for the final integration over z.
    part1_integrand_z = 2 * z * (1 - z/4)**6
    part2_integrand_z = (z/2) * (1 - z/4)

    # Integrate both parts with respect to z from 0 to 4
    flux_part1 = sympy.integrate(part1_integrand_z, (z, 0, 4))
    flux_part2 = sympy.integrate(part2_integrand_z, (z, 0, 4))

    # The total flux is 2 * (flux from part 1 + flux from part 2)
    flux_one_side = flux_part1 + flux_part2
    total_flux = 2 * flux_one_side
    
    # --- Output the Explanation and Result ---
    print("The total energy flow is twice the flow through one yellow side due to symmetry.")
    print("The flow through one side is calculated by the surface integral:")
    print("Flux_one_side = \u222B[z=0 to 4] \u222B[y=-(1-z/4) to (1-z/4)] (3zy\u00B2(1-z/4)\u00B3 + z/4) dy dz\n")
    print("This integral can be solved in two main steps.")
    print(f"1. Integrate \u222B 2z(1-z/4)\u2076 dz from 0 to 4. Result = {flux_part1.p}/{flux_part1.q}")
    print(f"2. Integrate \u222B (z/2)(1-z/4) dz from 0 to 4. Result = {flux_part2.p}/{flux_part2.q}\n")
    print("The flow through one side is the sum of these two parts.")
    print(f"Flow_one_side = {flux_part1.p}/{flux_part1.q} + {flux_part2.p}/{flux_part2.q} = {flux_one_side.p}/{flux_one_side.q}\n")
    print("The total flow through both yellow sides is twice this amount.")
    print("The final equation is:")
    print(f"Total Flow = 2 * ( {flux_part1.p}/{flux_part1.q} + {flux_part2.p}/{flux_part2.q} ) = {total_flux.p}/{total_flux.q}")


if __name__ == '__main__':
    calculate_energy_flow()
