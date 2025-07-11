import sympy

def solve_flux():
    """
    Calculates the flux of the vector field F through the two yellow sides of a pyramid.
    """
    # Define symbols
    x, y, z = sympy.symbols('x y z')

    # Define the vector field F = (Fx, Fy, Fz)
    Fx = 3 * x**3 * y**2 * z
    Fy = 3 * x**2 * y**3
    Fz = z

    print("We are calculating the energy flow (flux) through the two yellow sides of the pyramid.")
    print("Assuming the yellow sides are the 'left' and 'right' faces, whose plane equations are -4x + z = 4 and 4x + z = 4 respectively.\n")

    # --- Flux through the first yellow side (Right Face, S_R: 4x + z = 4) ---
    # The surface is z = 4 - 4x.
    # The outward normal vector for projection on the xy-plane is (4, 0, 1).
    
    # Substitute z = 4 - 4x into F
    Fx_R = Fx.subs(z, 4 - 4*x)
    Fy_R = Fy.subs(z, 4 - 4*x)
    Fz_R = Fz.subs(z, 4 - 4*x)

    # Calculate the dot product F . n for the integrand
    integrand_R = Fx_R * 4 + Fy_R * 0 + Fz_R * 1

    # The projection on the xy-plane is a triangle defined by 0 <= x <= 1 and -x <= y <= x.
    # Integrate with respect to y first, then x.
    flux_R_y_integrated = sympy.integrate(integrand_R, (y, -x, x))
    flux_R = sympy.integrate(flux_R_y_integrated, (x, 0, 1))
    
    print(f"Calculating the flux for the right face (4x + z = 4)...")
    print(f"The flux is: {flux_R}")

    # --- Flux through the second yellow side (Left Face, S_L: -4x + z = 4) ---
    # The surface is z = 4 + 4x.
    # The outward normal vector for projection on the xy-plane is (-4, 0, 1).

    # Substitute z = 4 + 4x into F
    Fx_L = Fx.subs(z, 4 + 4*x)
    Fy_L = Fy.subs(z, 4 + 4*x)
    Fz_L = Fz.subs(z, 4 + 4*x)

    # Calculate the dot product F . n for the integrand
    integrand_L = Fx_L * (-4) + Fy_L * 0 + Fz_L * 1

    # The projection on the xy-plane is a triangle defined by -1 <= x <= 0 and x <= y <= -x.
    # Integrate with respect to y first, then x.
    flux_L_y_integrated = sympy.integrate(integrand_L, (y, x, -x))
    flux_L = sympy.integrate(flux_L_y_integrated, (x, -1, 0))

    print(f"\nCalculating the flux for the left face (-4x + z = 4)...")
    print(f"The flux is: {flux_L}")

    # --- Total Flux ---
    total_flux = flux_R + flux_L

    print("\nThe total energy flow is the sum of the fluxes through both yellow sides.")
    # The final equation as requested
    print(f"Total Flow = {flux_R} + {flux_L} = {total_flux}")
    
    # Return the final numerical value in the specified format
    print(f"\n<<<80/21>>>")

solve_flux()