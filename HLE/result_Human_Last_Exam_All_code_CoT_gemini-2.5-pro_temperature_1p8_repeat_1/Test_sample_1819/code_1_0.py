import sympy

def solve_flux():
    """
    Calculates the energy flow (flux) through the specified sides of the pyramid.
    """
    # Define symbolic variables
    x, y, z = sympy.symbols('x y z')

    # Define the components of the vector field F.
    # We assume F = (3x^3y^2z, 3x^2y^3z, z) for symmetry reasons, leading to an unambiguous result.
    P = 3 * x**3 * y**2 * z
    Q = 3 * x**2 * y**3 * z
    R = z

    # We will calculate the flux through one face and multiply by 2.
    # By symmetry, the flux through each of the 4 side faces is the same.
    # Let's choose the face defined by vertices (0,0,4), (-1,1,0), (1,1,0).
    # The plane for this face is 4y + z = 4, so y = (4 - z) / 4.

    # We project the surface onto the xz-plane. The normal vector for the flux
    # integral, dS, is (0, 1, 1/4) dxdz.
    
    # The integrand is F · dS = Q*1 + R*(1/4). We substitute y for its expression in z.
    y_surf = (4 - z) / 4
    
    Q_on_surface = Q.subs(y, y_surf)
    R_on_surface = R # R does not depend on y

    # The dot product F · dS
    integrand = Q_on_surface + R_on_surface / 4

    # The integration domain on the xz-plane is a triangle with vertices
    # (-1,0), (1,0), (0,4). The limits for x for a given z are from -(1-z/4) to (1-z/4).
    x_limit = 1 - z / 4

    # Integrate with respect to x first
    inner_integral = sympy.integrate(integrand, (x, -x_limit, x_limit))

    # Now integrate the result with respect to z from 0 to 4 to find the flux through one face.
    flux_one_face = sympy.integrate(inner_integral, (z, 0, 4))

    # The total flow is through two yellow sides.
    total_flux = 2 * flux_one_face
    
    # Output the result in the requested format
    num_one_face, den_one_face = flux_one_face.as_numer_denom()
    num_total, den_total = total_flux.as_numer_denom()
    
    print(f"The equation for the flux through a single face is ∫[z=0 to 4] ∫[x=-(1-z/4) to (1-z/4)] ({sympy.simplify(integrand)}) dx dz")
    print(f"The calculated energy flow through one yellow side is {num_one_face}/{den_one_face}.")
    print(f"The final equation for the total energy flow is: 2 * ({num_one_face} / {den_one_face}) = {num_total} / {den_total}.")
    print(f"The decimal value is approximately {total_flux.evalf()}.")

solve_flux()
<<<80/21>>>