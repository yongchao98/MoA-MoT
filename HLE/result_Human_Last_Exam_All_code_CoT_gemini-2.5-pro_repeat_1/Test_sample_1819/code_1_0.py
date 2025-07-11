import sympy

def solve_flux():
    """
    Calculates the energy flow through the yellow sides of a square pyramid.
    """
    # Define symbolic variables
    x, y, z = sympy.symbols('x y z')

    # Define the vector field F
    F_x = 3 * x**3 * y**2 * z
    F_y = 3 * x**2 * y**3
    F_z = z
    F = sympy.Matrix([F_x, F_y, F_z])

    print("The energy flow vector field is F = (3x^3*y^2*z, 3x^2*y^3, z).")
    print("The pyramid has base vertices at (+-1, +-1, 0) and apex at (0, 0, 4).")
    print("The yellow sides are chosen to be the left and right faces.")

    # --- Step 1: Flux through the right yellow face (S_R) ---
    print("\n--- Calculating Flux through the Right Yellow Face (S_R) ---")
    # Plane equation: 4x + z = 4  => x = 1 - z/4
    x_R = 1 - z/4
    print(f"The plane for S_R is 4x + z = 4, or x = {x_R}.")

    # Substitute x into F
    F_on_SR = F.subs(x, x_R)

    # Normal vector dS for projection on yz-plane
    # For a surface x=h(y,z), dS = (1, -h_y, -h_z)dydz
    # h(y,z) = 1 - z/4. h_y=0, h_z=-1/4.
    # dS = (1, 0, 1/4)dydz. This is an outward normal (x>0).
    dS_R = sympy.Matrix([1, 0, sympy.Rational(1, 4)])
    print(f"The outward normal vector differential dS_R is ({dS_R[0]}, {dS_R[1]}, {dS_R[2]})dydz.")

    # Dot product F . dS
    dot_prod_R = F_on_SR.dot(dS_R)
    print(f"The dot product F . dS_R is: {sympy.simplify(dot_prod_R)}")

    # Integration limits
    # z from 0 to 4
    # y from -(1-z/4) to (1-z/4)
    y_lim_lower = -(1 - z/4)
    y_lim_upper = 1 - z/4
    
    # Integrate with respect to y first
    integral_y_R = sympy.integrate(dot_prod_R, (y, y_lim_lower, y_lim_upper))
    
    # Integrate with respect to z
    flux_R = sympy.integrate(integral_y_R, (z, 0, 4))
    print(f"The flux through the right yellow face is Phi_R = integral from z=0 to 4 of (integral from y=-(1-z/4) to (1-z/4) of (F . dS_R) dy) dz = {flux_R}")

    # --- Step 2: Flux through the left yellow face (S_L) ---
    print("\n--- Calculating Flux through the Left Yellow Face (S_L) ---")
    # Plane equation: -4x + z = 4  => x = (z - 4)/4
    x_L = (z - 4) / 4
    print(f"The plane for S_L is -4x + z = 4, or x = {x_L}.")

    # Substitute x into F
    F_on_SL = F.subs(x, x_L)
    
    # Normal vector dS for projection on yz-plane
    # For x=h(y,z), the normal is (1, -h_y, -h_z).
    # h(y,z) = (z-4)/4. h_y=0, h_z=1/4. Normal is (1,0,-1/4).
    # This points right (x>0). We need left-pointing normal for this face (x<0).
    # So we take the negative: dS_L = (-1, 0, 1/4)dydz.
    dS_L = sympy.Matrix([-1, 0, sympy.Rational(1, 4)])
    print(f"The outward normal vector differential dS_L is ({dS_L[0]}, {dS_L[1]}, {dS_L[2]})dydz.")

    # Dot product F . dS
    dot_prod_L = F_on_SL.dot(dS_L)
    print(f"The dot product F . dS_L is: {sympy.simplify(dot_prod_L)}")

    # Integrate with respect to y first
    integral_y_L = sympy.integrate(dot_prod_L, (y, y_lim_lower, y_lim_upper))
    
    # Integrate with respect to z
    flux_L = sympy.integrate(integral_y_L, (z, 0, 4))
    print(f"The flux through the left yellow face is Phi_L = {flux_L}")

    # --- Step 3: Total flux ---
    total_flux = flux_R + flux_L
    print("\n--- Total Flux through Yellow Sides ---")
    print(f"The total energy flow is the sum of the fluxes through the two yellow sides.")
    print(f"Total Flux = Phi_R + Phi_L = {flux_R} + {flux_L} = {total_flux}")
    
solve_flux()
print(f"\nFinal Answer: The total energy flow through the yellow sides is 80/21.")
<<<80/21>>>