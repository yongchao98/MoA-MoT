import sympy

def calculate_pyramid_flux():
    """
    Calculates the flux of a vector field F through the two yellow sides of a pyramid.
    """
    x, y, z = sympy.symbols('x y z')

    # Define the vector field F
    # F = (3x^3*y^2*z, 3x^2*y^3, z)
    
    # --- Step 1: Flux through the Front Yellow Face (S1) ---
    # This face is in the region y > 0.
    # The plane equation for this face is z = 4 - 4y.
    # The normal vector element dS is (0, 4, 1) * dA, pointing outwards.
    # The dot product F . dS = (3x^3*y^2*z, 3x^2*y^3, z) . (0, 4, 1)
    # = 12*x^2*y^3 + z.
    # Substitute z = 4 - 4y.
    integrand1_expr = 12 * x**2 * y**3 + (4 - 4 * y)
    
    # The projection on the xy-plane is a triangle with vertices (0,0), (-1,1), (1,1).
    # Integration bounds: 0 <= y <= 1, and -y <= x <= y.
    
    # Integrate with respect to x first
    inner_integral1 = sympy.integrate(integrand1_expr, (x, -y, y))
    
    # Integrate with respect to y
    flux1 = sympy.integrate(inner_integral1, (y, 0, 1))
    
    # --- Step 2: Flux through the Back Yellow Face (S2) ---
    # This face is in the region y < 0. By symmetry analysis, the flux is identical to the front face.
    # Let's verify by calculation.
    # The plane equation is z = 4 + 4y.
    # The normal vector element dS is (0, -4, 1) * dA, pointing outwards.
    # F . dS = (3x^3*y^2*z, 3x^2*y^3, z) . (0, -4, 1)
    # = -12*x^2*y^3 + z.
    # Substitute z = 4 + 4y.
    integrand2_expr = -12 * x**2 * y**3 + (4 + 4 * y)
    
    # The projection on the xy-plane is a triangle with vertices (0,0), (-1,-1), (1,-1).
    # Integration bounds: -1 <= y <= 0, and y <= x <= -y.

    # Integrate with respect to x first
    inner_integral2 = sympy.integrate(integrand2_expr, (x, y, -y))
    
    # Integrate with respect to y
    flux2 = sympy.integrate(inner_integral2, (y, -1, 0))

    # --- Step 3: Total Flux ---
    total_flux = flux1 + flux2
    
    print("Plan:")
    print("1. Calculate flux through the front yellow face (S1).")
    print("2. Calculate flux through the back yellow face (S2).")
    print("3. Sum the results for the total energy flow.\n")
    
    print("Calculation for the front face (S1):")
    print(f"Integral: ∫(from 0 to 1) ∫(from -y to y) ({integrand1_expr}) dx dy")
    print(f"Flux through S1 = {flux1}\n")

    print("Calculation for the back face (S2):")
    print(f"Integral: ∫(from -1 to 0) ∫(from y to -y) ({integrand2_expr}) dx dy")
    print(f"Flux through S2 = {flux2}\n")

    print("Total energy flow through the yellow sides:")
    print(f"{flux1} + {flux2} = {total_flux}")

if __name__ == '__main__':
    calculate_pyramid_flux()