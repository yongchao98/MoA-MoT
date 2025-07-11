import sympy as sp

def solve_flux():
    """
    This function calculates the energy flow through the yellow sides of the pyramid.
    """
    # Define the symbols for our coordinate system
    x, y, z = sp.symbols('x y z')

    # Define the vector field F based on the problem description
    # F = (3x^3*y^2*z, 3x^2*y^3*z, z)
    F = sp.Matrix([3 * x**3 * y**2 * z, 3 * x**2 * y**3 * z, z])

    # We calculate the flux for one of the four identical slanted faces.
    # Let's choose the face in the back (positive y), which lies on the plane 4y + z = 4.
    # We parameterize this surface by projecting it onto the xz-plane.
    # From the plane equation, we get y = 1 - z/4.
    # The parameterization r(x, z) is:
    r = sp.Matrix([x, 1 - z/4, z])

    # The domain for the parameters (x, z) is the triangle with vertices
    # (-1, 0), (1, 0), and (0, 4).
    # We can exploit the symmetry of the setup. We'll integrate x from 0 to 1
    # and z from 0 up to the line defined by the vertices (1,0) and (0,4),
    # which is z = 4 - 4x. We then multiply the result by 2 to account for the
    # symmetric half (from x = -1 to 0).

    # Compute the partial derivatives of the parameterization
    r_diff_x = r.diff(x)
    r_diff_z = r.diff(z)

    # The cross product gives a normal vector to the surface.
    # We need the outward-pointing normal. The gradient of the plane g=4y+z-4 is (0,4,1),
    # which points outward. The cross product r_diff_x.cross(r_diff_z) is (0, -1, -1/4),
    # which points inward. So we take its negative.
    outward_normal = -r_diff_x.cross(r_diff_z)

    # Substitute the parameterization for y into the vector field F
    F_on_surface = F.subs(y, r[1])

    # The integrand is the dot product of F on the surface and the outward normal
    integrand = F_on_surface.dot(outward_normal)

    # Integrate with respect to z first, over the range [0, 4 - 4x]
    integral_over_z = sp.integrate(integrand, (z, 0, 4 - 4*x))

    # Integrate the result with respect to x over [0, 1]
    flux_for_half_face = sp.integrate(integral_over_z, (x, 0, 1))

    # The total flux for one face is twice the result due to symmetry
    flux_one_side = 2 * flux_for_half_face
    
    # As reasoned in the plan, all four slanted faces have an identical flux value.
    # Since two sides are painted yellow, the total flow is twice the flow through one side.
    total_flux = 2 * flux_one_side

    # Print the step-by-step results
    print(f"The energy flow (flux) through a single yellow side is calculated by a surface integral.")
    print(f"The result of this integral for one side is: {flux_one_side}")
    print(f"Since there are two identical yellow sides, we sum their contributions.")
    print(f"Total Flow = {flux_one_side} + {flux_one_side} = {total_flux}")
    print(f"The numerical value of the total flow is approximately: {total_flux.evalf()}")

solve_flux()
<<<80/21>>>