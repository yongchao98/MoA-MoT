import sympy

def solve_flux():
    """
    Calculates the flux of a vector field F through the yellow sides of a pyramid.
    """
    # Define the symbols
    x, y, z = sympy.symbols('x y z')

    # Define the vector field F = (Fx, Fy, Fz)
    Fx = 3 * x**3 * y**2 * z
    Fy = 3 * x**2 * y**3
    Fz = z

    # --- Step 1: Calculate the flux through the right yellow side ---
    # This side lies on the plane 4x + z = 4, or x = 1 - z/4.
    # The outward normal vector is (4, 0, 1).
    # The differential surface element dS for projection onto the yz-plane is (1, 0, 1/4) dy dz.
    # The dot product F . dS is (Fx * 1) + (Fy * 0) + (Fz * 1/4).
    # We substitute the plane equation x = 1 - z/4 into this dot product.
    
    integrand_right = Fx.subs(x, 1 - z/4) + Fz / 4

    # The projection of this triangular face onto the yz-plane is a triangle
    # with vertices (y,z) at (1,0), (-1,0), and (0,4).
    # The integration limits are z from 0 to 4, and y from -(1-z/4) to (1-z/4).

    # Integrate with respect to y first
    inner_integral = sympy.integrate(integrand_right, (y, -(1 - z/4), 1 - z/4))

    # Now integrate the result with respect to z to get the flux for one side
    flux_one_side = sympy.integrate(inner_integral, (z, 0, 4))

    # --- Step 2: Calculate the total flux ---
    # Due to symmetry, the flux through the left yellow side is identical.
    # Total flux is twice the flux of one side.
    total_flux = 2 * flux_one_side

    # --- Step 3: Print the results ---
    print("The energy flow through the yellow sides is calculated by summing the flux through each side.")
    print("The vector field is F = (3x^3*y^2*z, 3x^2*y^3, z).")
    print("The pyramid side in the x > 0 region is on the plane 4x + z = 4.")
    print(f"The flux through one yellow side is calculated by the integral:")
    print(f"Flux_one_side = Integral from z=0 to 4 of (Integral from y=-(1-z/4) to (1-z/4) of ({integrand_right}) dy) dz")
    
    # Extract numerator and denominator for clean printing
    num, den = sympy.fraction(total_flux)
    
    print("\n--- Calculation Result ---")
    print(f"Flux through the first yellow side = {flux_one_side}")
    print(f"Flux through the second yellow side = {flux_one_side}")
    print("Total energy flow = Flux_side_1 + Flux_side_2")
    print(f"Total energy flow = {flux_one_side} + {flux_one_side} = {total_flux}")
    
solve_flux()
<<<80/21>>>