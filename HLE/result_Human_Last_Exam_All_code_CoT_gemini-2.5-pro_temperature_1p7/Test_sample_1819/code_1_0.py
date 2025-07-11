import sympy as sp

def solve_flux():
    """
    Calculates the flux of a vector field F through two opposite sides of a square pyramid.
    """
    # Define symbols
    x, y, z = sp.symbols('x y z')

    # Define the vector field F
    Fx = 3 * x**3 * y**2 * z
    Fy = 3 * x**2 * y**3
    Fz = z
    F = sp.Matrix([Fx, Fy, Fz])

    print("Step 1: Define the vector field and the surface of one yellow side.")
    print("We assume the yellow sides are the opposite faces in the x>0 and x<0 regions.")
    print(f"Vector Field F = ({Fx}, {Fy}, {Fz})")
    
    # The face in the x > 0 region is on the plane 4x + z = 4.
    # We can parameterize this surface with x = 1 - z/4.
    # The parameterization vector is r(y, z) = [1 - z/4, y, z].
    
    # Create the parameterization vector
    r = sp.Matrix([1 - z/4, y, z])

    # Find the tangent vectors with respect to y and z
    r_y = sp.diff(r, y)
    r_z = sp.diff(r, z)

    # Compute the normal vector to the surface.
    # The cross product r_y x r_z gives an outward-pointing normal vector.
    N = r_y.cross(r_z)
    
    # Substitute the parameterization for x into the vector field F
    F_on_surface = F.subs(x, 1 - z/4)

    # Calculate the dot product F · N for the integrand
    integrand = F_on_surface.dot(N)
    
    print("\nStep 2: Set up the surface integral.")
    print(f"The integrand (F · dS) is: {integrand}")
    
    # The integration limits are z from 0 to 4, and for each z, y goes from -(1 - z/4) to (1 - z/4).
    # Integrate with respect to y first
    inner_integral = sp.integrate(integrand, (y, -(1 - z/4), (1 - z/4)))

    # Integrate the result with respect to z to find the flux
    flux_one_side = sp.integrate(inner_integral, (z, 0, 4))
    
    print("\nStep 3: Calculate the flux for one yellow side.")
    print(f"The result of the double integral for one face is: {flux_one_side}")
    
    # The flux through the opposite face is identical.
    total_flux = 2 * flux_one_side
    
    print("\nStep 4: Calculate the total flux for the two yellow sides.")
    print(f"The total energy flow is the sum of the fluxes through both yellow sides.")
    print(f"Total Flow = {flux_one_side} + {flux_one_side} = {total_flux}")

if __name__ == '__main__':
    solve_flux()
    # The final answer in the required format
    print("\n<<<80/21>>>")
