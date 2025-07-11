import sympy

def solve_flux():
    """
    Calculates the flux of a vector field F through the two yellow sides of a pyramid.
    """
    # Define symbolic variables
    x, y, z = sympy.symbols('x y z')

    # Define the vector field F = (P, Q, R)
    F_P = 3 * x**3 * y**2 * z
    F_Q = 3 * x**2 * y**3
    F_R = z
    F = sympy.Matrix([F_P, F_Q, F_R])

    # --- Calculation for the first yellow side (S1, on x > 0) ---

    # S1 is on the plane z = 4 - 4x
    z1_expr = 4 - 4*x
    
    # Parameterize S1 by projecting on the xy-plane.
    # The normal vector for a surface z=f(x,y) is (-dz/dx, -dz/dy, 1).
    # This must be the outward-pointing normal.
    # dz/dx = -4, dz/dy = 0. So normal is (4, 0, 1). This points outwards (+x direction).
    n1 = sympy.Matrix([4, 0, 1])

    # F on the surface S1
    F_on_S1 = F.subs(z, z1_expr)

    # Dot product F . n
    dot_prod1 = F_on_S1.dot(n1)

    # The projection of S1 on the xy-plane is a triangle with vertices (0,0), (1,-1), (1,1).
    # The integration limits are 0 <= x <= 1, -x <= y <= x.
    # Integrate to find the flux
    flux1 = sympy.integrate(dot_prod1, (y, -x, x), (x, 0, 1))

    # --- Calculation for the second yellow side (S2, on x < 0) ---

    # S2 is on the plane z = 4 + 4x
    z2_expr = 4 + 4*x

    # The outward-pointing normal is (-dz/dx, -dz/dy, 1) = (-4, 0, 1).
    n2 = sympy.Matrix([-4, 0, 1])

    # F on the surface S2
    F_on_S2 = F.subs(z, z2_expr)
    
    # Dot product F . n
    dot_prod2 = F_on_S2.dot(n2)
    
    # The projection of S2 on the xy-plane is a triangle with vertices (0,0), (-1,1), (-1,-1).
    # The integration limits are -1 <= x <= 0, x <= y <= -x.
    flux2 = sympy.integrate(dot_prod2, (y, x, -x), (x, -1, 0))

    # --- Total Flux ---
    total_flux = flux1 + flux2
    
    print("The vector field is F = (3*x^3*y^2*z, 3*x^2*y^3, z)")
    print("\nThe first yellow side's flux calculation:")
    print(f"Flux_1 = {flux1}")
    
    print("\nThe second yellow side's flux calculation:")
    print(f"Flux_2 = {flux2}")
    
    print("\nThe total energy flow is the sum of the fluxes through the two yellow sides.")
    # The output format requires printing the equation with numbers
    print(f"Total Flux = {flux1} + {flux2} = {total_flux}")
    
solve_flux()