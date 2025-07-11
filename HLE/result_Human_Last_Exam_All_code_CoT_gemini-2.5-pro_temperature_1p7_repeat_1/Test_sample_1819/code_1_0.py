import sympy as sp

def solve_flux():
    """
    Calculates the flux of the vector field F through the two yellow sides of the pyramid.
    """
    # Define symbols
    x, y, z = sp.symbols('x y z')

    # Define the vector field F = (P, Q, R)
    P = 3 * x**3 * y**2 * z
    Q = 3 * x**2 * y**3
    R = z

    # --- Flux through the Front Face (S1) ---
    # The surface S1 is defined by the plane z = 4 - 4x.
    # We use the formula for flux over a surface z=g(x,y):
    # Flux = integral(-P*gx - Q*gy + R) dA, where gx and gy are partial derivatives of g.
    # Here g(x,y) = 4 - 4x, so gx = -4 and gy = 0.
    # The normal vector points outwards.
    
    # Integrand for S1, substituting z = 4 - 4x
    integrand1 = (-P * (-4) - Q * 0 + R).subs(z, 4 - 4*x)
    
    # The projection of S1 on the xy-plane is a triangle with vertices (0,0), (1,-1), (1,1).
    # The integration limits are x from 0 to 1, and y from -x to x.
    flux1 = sp.integrate(integrand1, (y, -x, x), (x, 0, 1))

    # --- Flux through the Back Face (S3) ---
    # The surface S3 is defined by the plane z = 4 + 4x.
    # g(x,y) = 4 + 4x, so gx = 4 and gy = 0.
    
    # Integrand for S3, substituting z = 4 + 4x
    integrand3 = (-P * 4 - Q * 0 + R).subs(z, 4 + 4*x)
    
    # The projection of S3 on the xy-plane is a triangle with vertices (0,0), (-1,-1), (-1,1).
    # The integration limits are x from -1 to 0, and y from x to -x.
    # Due to symmetry, the flux through S3 is the same as through S1.
    flux3 = sp.integrate(integrand3, (y, x, -x), (x, -1, 0))

    # --- Total Flux ---
    total_flux = flux1 + flux3

    print("We calculate the energy flow, which is the flux of the vector field F through the two yellow sides.")
    print("Assuming the 'front' and 'back' faces of the pyramid are yellow:")
    print(f"The flux through the front face is calculated as: {flux1}")
    print(f"The flux through the back face is calculated as: {flux3}")
    print("\nThe total energy flow is the sum of the fluxes through these two faces.")
    print(f"Total Flow = {flux1} + {flux3} = {total_flux}")

solve_flux()
<<<80/21>>>