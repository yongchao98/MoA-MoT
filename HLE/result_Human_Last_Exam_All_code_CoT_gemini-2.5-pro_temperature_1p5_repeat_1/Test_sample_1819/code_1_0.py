import sympy

def solve_flux():
    """
    Calculates the flux of a vector field F through the two yellow sides of a pyramid.
    """
    # Define symbols
    x, y, z = sympy.symbols('x y z')

    # Define the vector field F
    F = sympy.Matrix([3 * x**3 * y**2 * z, 3 * x**2 * y**3, z])

    # --- Yellow Face 1 (y > 0) ---
    # Plane equation: z = 4 - 4y
    # Parameterization over xy-plane. Projection is a triangle with vertices (0,0), (-1,1), (1,1).
    # Integration bounds: y from 0 to 1, x from -y to y.
    
    # Substitute z into F for this surface
    z_s1 = 4 - 4*y
    F_s1 = F.subs(z, z_s1)
    
    # Outward normal vector dS = (0, 4, 1) dx dy
    n1 = sympy.Matrix([0, 4, 1])
    
    # Integrand is the dot product F . dS
    integrand1 = F_s1.dot(n1)
    
    # Integrate over the domain: first x, then y
    flux1 = sympy.integrate(integrand1, (x, -y, y), (y, 0, 1))

    # --- Yellow Face 2 (y < 0) ---
    # Plane equation: z = 4 + 4y
    # Parameterization over xy-plane. Projection is a triangle with vertices (0,0), (1,-1), (-1,-1).
    # Integration bounds: y from -1 to 0, x from y to -y.

    # Substitute z into F for this surface
    z_s2 = 4 + 4*y
    F_s2 = F.subs(z, z_s2)
    
    # Outward normal vector dS = (0, -4, 1) dx dy
    n2 = sympy.Matrix([0, -4, 1])
    
    # Integrand is the dot product F . dS
    integrand2 = F_s2.dot(n2)
    
    # Integrate over the domain: first x, then y
    flux2 = sympy.integrate(integrand2, (x, y, -y), (y, -1, 0))

    # Total flux is the sum of the two
    total_flux = flux1 + flux2

    # Print the results step-by-step
    print(f"The flux through the first yellow face (y>0) is calculated to be: {flux1}")
    print(f"The flux through the second yellow face (y<0) is calculated to be: {flux2}")
    print("\nThe total energy flow through the yellow sides is the sum of these two values:")
    print(f"Total Flow = {flux1} + {flux2} = {total_flux}")

solve_flux()
<<<104/21>>>