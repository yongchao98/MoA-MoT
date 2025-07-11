import sympy

def solve_flux():
    """
    Calculates the flux of the vector field F through the yellow sides of a pyramid.
    """
    # Define symbolic variables
    x, y, z = sympy.symbols('x y z')

    # Define the vector field F = (P, Q, R)
    P = 3 * x**3 * y**2 * z
    Q = 3 * x**2 * y**3
    R = z
    
    # The pyramid's side faces are defined by the planes connecting the base to the apex (0,0,4).
    # The base vertices are at (±1, ±1, 0).
    # The projection of the side faces onto a coordinate plane is a triangle. For a given height z,
    # the cross-section is a square with side length 2*(1-z/4).
    # Let L(z) = 1 - z/4. At height z, x and y are within [-L(z), L(z)].

    # We assume the yellow sides are the front and back faces.

    # --- Flux through the Front Face (y > 0) ---
    # The front face lies on the plane y = 1 - z/4.
    # The outward normal vector component for the flux integral is (0, 1, 1/4).
    # The flux integral is ∫∫ (Q + R/4) dx dz over the projection on the xz-plane.
    y_front = 1 - z/4
    integrand_front = Q.subs(y, y_front) + R / 4
    L = 1 - z/4
    
    # Integrate with respect to x first, then z
    flux_F_inner = sympy.integrate(integrand_front, (x, -L, L))
    flux_F = sympy.integrate(flux_F_inner, (z, 0, 4))

    # --- Flux through the Back Face (y < 0) ---
    # The back face lies on the plane y = -(1 - z/4).
    # The outward normal vector component is (0, -1, 1/4).
    # The flux integral is ∫∫ (-Q + R/4) dx dz.
    y_back = -(1 - z/4)
    integrand_back = -Q.subs(y, y_back) + R / 4
    
    # Integrate with respect to x first, then z
    flux_B_inner = sympy.integrate(integrand_back, (x, -L, L))
    flux_B = sympy.integrate(flux_B_inner, (z, 0, 4))

    # Total flux is the sum of the fluxes through the two yellow faces
    total_flux = flux_F + flux_B

    print("Assuming the yellow sides are the front and back faces of the pyramid.")
    print("The total energy flow is the sum of the fluxes through these two faces.")
    print(f"Flux through front face = {flux_F}")
    print(f"Flux through back face = {flux_B}")
    print(f"Total flow through yellow sides = {flux_F} + {flux_B} = {total_flux}")

solve_flux()