import sympy

def solve_flux():
    """
    Calculates the energy flow (flux) through the two specified yellow sides of the pyramid.
    """
    # Define symbolic variables
    x, y, z = sympy.symbols('x y z')

    # Define the vector field F
    F = sympy.Matrix([3 * x**3 * y**2 * z, 3 * x**2 * y**3, z])

    # We select the front face (in plane 4y + z = 4) and back face (-4y + z = 4) as the yellow sides.
    # --- Calculation for the Front Face ---
    # The surface normal is dS = (0, 1, 1/4) dx dz. On this face, y = 1 - z/4.
    y_front = 1 - z/4
    dS_front = sympy.Matrix([0, 1, sympy.Rational(1, 4)])

    # Substitute y into F and compute the dot product for the integrand
    F_on_front = F.subs(y, y_front)
    integrand_front = F_on_front.dot(dS_front)

    # Define integration limits for the projection on the xz-plane
    # z from 0 to 4, x from -(1-z/4) to (1-z/4)
    x_lim = 1 - z/4
    
    # Calculate the flux for the front face
    flux_front = sympy.integrate(integrand_front, (x, -x_lim, x_lim), (z, 0, 4))

    # --- Calculation for the Back Face ---
    # The integrand and integration limits for the back face turn out to be identical.
    # Therefore, the flux through the back face is the same.
    flux_back = flux_front
    
    # The total flux is the sum of the fluxes through the two yellow sides.
    total_flux = flux_front + flux_back

    # Print the final equation
    print(f"The energy flow through one yellow side is {flux_front}.")
    print(f"The energy flow through the other yellow side is {flux_back}.")
    print(f"The total energy flow through the yellow sides is {flux_front} + {flux_back} = {total_flux}")

solve_flux()