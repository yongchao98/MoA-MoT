import sympy as sp

def solve_flux():
    """
    This function calculates the total energy flow (flux) through the two yellow sides of a pyramid.
    """
    # Define symbols for our variables
    x, y, z = sp.symbols('x y z')

    # Define the vector field F
    F_vector = sp.Matrix([3 * x**3 * y**2 * z, 3 * x**2 * y**3, z])

    # --- Plan ---
    # We will calculate the flux for the two opposite sides which we define as "yellow".
    # Let's choose the faces in the positive and negative x directions.
    # Yellow Side 1 (S1): The face on the plane 4x + z = 4.
    # Yellow Side 2 (S2): The face on the plane -4x + z = 4.
    
    # --- Calculation for Yellow Side 1 (S1: 4x + z = 4) ---
    
    # To calculate the flux integral ∫∫_S F · dS, we can project the surface onto the yz-plane.
    # The surface is given by x = g(y, z) = 1 - z/4.
    # The differential surface element dS is given by n_vec * dy * dz, where n_vec is (1, -∂g/∂y, -∂g/∂z).
    # ∂g/∂y = 0 and ∂g/∂z = -1/4. So, n_vec = (1, 0, 1/4). This vector points in the +x direction, which is outward.
    
    n_vec1 = sp.Matrix([1, 0, sp.Rational(1, 4)])
    
    # Substitute x on the surface into the vector field F
    F_on_S1 = F_vector.subs(x, 1 - z/4)
    
    # The integrand is the dot product of F on the surface and the normal vector.
    integrand1 = F_on_S1.dot(n_vec1)
    
    # The surface projects onto a triangle in the yz-plane with vertices (0,4), (1,0), and (-1,0).
    # The limits of integration are z from 0 to 4, and for each z, y from -(1 - z/4) to (1 - z/4).
    y_limit = 1 - z/4
    
    # Perform the double integration
    flux1 = sp.integrate(integrand1, (y, -y_limit, y_limit), (z, 0, 4))
    
    # --- Calculation for Yellow Side 2 (S2: -4x + z = 4) ---
    
    # The surface is given by x = g(y, z) = z/4 - 1.
    # The normal vector n_vec = (1, -∂g/∂y, -∂g/∂z) = (1, 0, -1/4), which points in the +x direction.
    # For this face, the outward direction is -x, so we use the opposite normal vector.
    n_vec2 = sp.Matrix([-1, 0, sp.Rational(1, 4)])
    
    # Substitute x on the surface into the vector field F
    F_on_S2 = F_vector.subs(x, z/4 - 1)
    
    # The integrand is the dot product.
    integrand2 = F_on_S2.dot(n_vec2)
    
    # The limits of integration are the same.
    flux2 = sp.integrate(integrand2, (y, -y_limit, y_limit), (z, 0, 4))
    
    # The total flux is the sum of the fluxes through the two yellow sides.
    total_flux = flux1 + flux2

    # Print the results in a clear format
    print("The total energy flow is the sum of the fluxes through the two yellow sides.")
    print(f"The flux through the first yellow side is: {flux1}")
    print(f"The flux through the second yellow side is: {flux2}")
    print("The final equation for the total flow is:")
    print(f"Total Flow = {flux1} + {flux2} = {total_flux}")

solve_flux()
<<<80/21>>>