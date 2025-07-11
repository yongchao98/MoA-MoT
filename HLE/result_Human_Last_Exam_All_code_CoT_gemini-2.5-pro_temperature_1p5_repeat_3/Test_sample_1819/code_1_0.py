import sympy

def solve_flux():
    """
    Calculates the flux of a vector field through the yellow sides of a pyramid.
    """
    # Define the symbolic variables
    x, y, z = sympy.symbols('x y z')

    # Define the vector field F = (P, Q, R)
    F_P = 3 * x**3 * y**2 * z
    F_Q = 3 * x**2 * y**3
    F_R = z
    F = sympy.Matrix([F_P, F_Q, F_R])

    # --- Calculation for one yellow face ---
    # The first yellow face is a triangle with vertices (1,-1,0), (1,1,0), (0,0,4).
    # The plane containing this face is 4x + z = 4, or x = 1 - z/4.
    
    # We parameterize the surface using y and z: r(y,z) = <1-z/4, y, z>
    # The outward normal vector for a surface x=g(y,z) is <1, -dg/dy, -dg/dz>.
    # Here, g(y,z) = 1-z/4, so dg/dy = 0 and dg/dz = -1/4.
    # The normal vector N is <1, 0, 1/4>.
    
    # Substitute the parameterization into F
    x_param = 1 - z/4
    F_on_surface = F.subs({x: x_param})

    # The normal vector for integration
    N = sympy.Matrix([1, 0, sympy.Rational(1, 4)])
    
    # The integrand is the dot product of F on the surface and the normal vector N.
    # Integrand = (3*(1-z/4)**3 * y**2 * z) * 1 + (z) * (1/4)
    integrand = F_on_surface.dot(N)
    
    term1, term2 = integrand.as_add_args()
    
    # The integration domain is the projection of the face onto the yz-plane.
    # This is a triangle with vertices (-1,0), (1,0), (0,4).
    # Limits for z are [0, 4].
    # Limits for y are [z/4 - 1, 1 - z/4].
    y_lim_lower = z/4 - 1
    y_lim_upper = 1 - z/4
    z_lim_lower = 0
    z_lim_upper = 4

    # Integrate the first term of the integrand
    flux_part1 = sympy.integrate(term1, (y, y_lim_lower, y_lim_upper), (z, z_lim_lower, z_lim_upper))
    
    # Integrate the second term of the integrand
    flux_part2 = sympy.integrate(term2, (y, y_lim_lower, y_lim_upper), (z, z_lim_lower, z_lim_upper))
    
    # The flux through one yellow face is the sum of the two parts
    flux_one_face = flux_part1 + flux_part2
    
    # Due to symmetry, the flux through the second yellow face is identical.
    # The total flux is twice the flux through one face.
    total_flux = 2 * flux_one_face

    # Print the results step-by-step
    print(f"The integrand is composed of two terms. We integrate them separately.")
    print(f"Integral of the first term: {flux_part1}")
    print(f"Integral of the second term: {flux_part2}")
    print(f"Flux through one yellow face = {flux_part1} + {flux_part2} = {flux_one_face}")
    print(f"Total energy flow through both yellow sides = 2 * ({flux_one_face}) = {total_flux}")

solve_flux()