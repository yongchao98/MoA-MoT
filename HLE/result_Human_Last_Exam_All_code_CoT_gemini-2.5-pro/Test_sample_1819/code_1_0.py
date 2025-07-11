import sympy
from sympy import Rational

def solve_flux():
    """
    Calculates the flux of a vector field F through the yellow sides of a pyramid.
    """
    # Define symbols for our calculations
    x, y, z = sympy.symbols('x y z')

    # Define the vector field F = (3x^3*y^2*z, 3x^2*y^3, z)
    F = sympy.Matrix([3 * x**3 * y**2 * z, 3 * x**2 * y**3, z])

    print("Calculating the energy flow (flux) through the yellow sides of the pyramid.")
    print("Vector Field F = (3x^3*y^2*z, 3x^2*y^3, z)")
    print("The pyramid has its apex at (0,0,4) and base vertices at (+-1, +-1, 0).")
    print("We assume the two yellow faces are the front and back faces (base edges parallel to the x-axis).\n")

    # --- Flux through the front face (S1) ---
    # The front face is a triangle with vertices (1,1,0), (-1,1,0), and (0,0,4).
    # Its plane equation is 4y + z = 4, so z = 4 - 4y.
    # We parameterize the surface as r(x,y) = (x, y, 4-4y).
    # The outward normal vector is dS = (0, 4, 1) dx dy.
    print("--- Calculating flux through the front yellow face (S1) ---")
    z_s1 = 4 - 4 * y
    F_on_s1 = F.subs(z, z_s1)
    normal_s1 = sympy.Matrix([0, 4, 1])
    
    # The dot product F . dS gives the integrand
    integrand1 = F_on_s1.dot(normal_s1)
    
    # The integration domain D1 is the projection of the face onto the xy-plane,
    # which is a triangle with vertices (0,0), (1,1), and (-1,1).
    # This corresponds to integration bounds: 0 <= y <= 1, -y <= x <= y.
    
    # Integrate with respect to x first, then y
    flux1_inner = sympy.integrate(integrand1, (x, -y, y))
    flux1 = sympy.integrate(flux1_inner, (y, 0, 1))
    
    flux1_val = Rational(flux1)
    print(f"The equation for the front face is z = 4 - 4y.")
    print(f"The integrand (F . n) is: {integrand1}")
    print(f"The flux through the front face is: {flux1_val}\n")

    # --- Flux through the back face (S3) ---
    # The back face has vertices (1,-1,0), (-1,-1,0), and (0,0,4).
    # Its plane equation is -4y + z = 4, so z = 4 + 4y.
    # Parameterization: r(x,y) = (x, y, 4+4y).
    # Outward normal vector dS = (0, -4, 1) dx dy.
    print("--- Calculating flux through the back yellow face (S3) ---")
    z_s3 = 4 + 4 * y
    F_on_s3 = F.subs(z, z_s3)
    normal_s3 = sympy.Matrix([0, -4, 1])

    # The dot product F . dS gives the integrand
    integrand3 = F_on_s3.dot(normal_s3)

    # The integration domain D3 is the projection onto the xy-plane:
    # a triangle with vertices (0,0), (1,-1), and (-1,-1).
    # Bounds: -1 <= y <= 0, y <= x <= -y.

    # Integrate with respect to x first, then y
    flux3_inner = sympy.integrate(integrand3, (x, y, -y))
    flux3 = sympy.integrate(flux3_inner, (y, -1, 0))

    flux3_val = Rational(flux3)
    print(f"The equation for the back face is z = 4 + 4y.")
    print(f"The integrand (F . n) is: {integrand3}")
    print(f"The flux through the back face is: {flux3_val}\n")
    
    # --- Total Flux ---
    total_flux = flux1 + flux3
    total_flux_val = Rational(total_flux)
    
    print("--- Total Energy Flow ---")
    print("The total energy flow through the yellow sides is the sum of the two fluxes.")
    print(f"Total energy flow = {flux1_val} + {flux3_val} = {total_flux_val}")

solve_flux()