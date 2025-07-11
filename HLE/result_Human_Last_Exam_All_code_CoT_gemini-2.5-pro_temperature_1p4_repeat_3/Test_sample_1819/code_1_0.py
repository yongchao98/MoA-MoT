import sympy

def solve_flux_problem():
    """
    This function calculates the flux of a vector field through two sides of a pyramid.
    """
    # Define the symbols
    x, y, z = sympy.symbols('x y z')

    # Define the vector field F
    F = sympy.Matrix([3*x**3*y**2*z, 3*x**2*y**3, z])

    # --- Pyramid Geometry & Plan ---
    # Apex: (0, 0, 4)
    # Base vertices: (±1, ±1, 0)
    # We will calculate the flux through the front and back faces, assuming they are yellow.

    print("Step 1: Calculate the flux through the front yellow face.")
    # The front face connects vertices (1, -1, 0), (-1, -1, 0) to the apex (0, 0, 4).
    # The equation of the plane for this face is z = 4y + 4.
    # An outward-pointing normal vector n is (0, -4, 1).
    z_front_eq = 4*y + 4
    n_front = sympy.Matrix([0, -4, 1])
    
    # Calculate the dot product F · n
    F_dot_n_front = F.dot(n_front)
    
    # Substitute z from the plane equation into the dot product
    integrand_front = F_dot_n_front.subs(z, z_front_eq)
    
    # The projection of the front face onto the xy-plane is a triangle
    # with vertices (0,0), (1,-1), (-1,-1).
    # The limits are: y from -1 to 0, and x from y to -y.
    flux_front = sympy.integrate(integrand_front, (x, y, -y), (y, -1, 0))
    print(f"The flux through the front face is: {flux_front}")
    
    print("\nStep 2: Calculate the flux through the back yellow face.")
    # The back face connects vertices (-1, 1, 0), (1, 1, 0) to the apex (0, 0, 4).
    # The equation of the plane for this face is z = 4 - 4y.
    # An outward-pointing normal vector n is (0, 4, 1).
    z_back_eq = 4 - 4*y
    n_back = sympy.Matrix([0, 4, 1])

    # Calculate the dot product F · n
    F_dot_n_back = F.dot(n_back)
    
    # Substitute z from the plane equation into the dot product
    integrand_back = F_dot_n_back.subs(z, z_back_eq)

    # The projection of the back face onto the xy-plane is a triangle
    # with vertices (0,0), (-1,1), (1,1).
    # The limits are: y from 0 to 1, and x from -y to y.
    flux_back = sympy.integrate(integrand_back, (x, -y, y), (y, 0, 1))
    print(f"The flux through the back face is: {flux_back}")

    print("\nStep 3: Calculate the total energy flow.")
    # The total flux is the sum of the fluxes through the two yellow faces.
    total_flux = flux_front + flux_back
    print("The total energy flow is the sum of the two fluxes.")
    print(f"Total Flux = {flux_front} + {flux_back} = {total_flux}")

solve_flux_problem()
<<<104/21>>>