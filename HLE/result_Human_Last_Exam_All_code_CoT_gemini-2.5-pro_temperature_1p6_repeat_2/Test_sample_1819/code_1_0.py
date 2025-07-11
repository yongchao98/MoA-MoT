import sympy

def solve_flux():
    """
    Calculates the flux of a vector field F through the two yellow sides of a pyramid.
    """
    # Define symbolic variables for parameterization
    u, v = sympy.symbols('u v')

    # --- Part 1: Flux through the Front Yellow Face (Y1) ---
    # The front face is a triangle with vertices: Apex(0,0,4), V1(1,1,0), V2(-1,1,0).
    # We can parameterize this surface as r(u,v) = (u*(2*v-1), u, 4-4*u) for u,v in [0,1].
    # This covers the triangle from the apex (u=0) to the base edge (u=1).

    # Components of the parameterization for the front face
    x1 = u * (2*v - 1)
    y1 = u
    z1 = 4 * (1 - u)

    # The vector field F = (3*x**3*y**2*z, 3*x**2*y**3, z). We evaluate it on the surface.
    # We only need the y and z components since the normal vector's x-component is 0.
    Fy1 = 3 * x1**2 * y1**3
    Fz1 = z1

    # The outward normal vector dS is calculated from the cross product of the partial
    # derivatives of the parameterization, dS = (r_v x r_u) du dv.
    # For the front face, this gives dS = (0, 8*u, 2*u) du dv.
    d_sy1 = 8 * u
    d_sz1 = 2 * u

    # The integrand is the dot product F . dS.
    # The x-component of F is multiplied by 0.
    integrand1 = Fy1 * d_sy1 + Fz1 * d_sz1
    integrand1_simplified = sympy.simplify(integrand1)

    # First, integrate with respect to v from 0 to 1
    integral_v1 = sympy.integrate(integrand1_simplified, (v, 0, 1))

    # Second, integrate the result with respect to u from 0 to 1 to get the flux
    flux1 = sympy.integrate(integral_v1, (u, 0, 1))

    # --- Part 2: Flux through the Back Yellow Face (Y2) ---
    # By symmetry, the calculation for the opposite back face yields the same flux.
    # The y-component of both F and the normal vector change sign, so their product remains the same.
    flux2 = flux1

    # --- Part 3: Total Flux ---
    # The total energy flow is the sum of the fluxes through the two yellow sides.
    total_flux = flux1 + flux2

    # --- Final Output ---
    print("To find the total energy flow through the yellow sides, we calculate the flux for each of the two faces and sum the results.")
    print("Assuming the front and back faces are yellow:")
    print(f"\nThe flux through the front yellow face is {flux1}.")
    print(f"The flux through the back yellow face is also {flux2}.")
    print("\nThe final equation for the total energy flow is:")
    print(f"{flux1} + {flux2} = {total_flux}")

solve_flux()