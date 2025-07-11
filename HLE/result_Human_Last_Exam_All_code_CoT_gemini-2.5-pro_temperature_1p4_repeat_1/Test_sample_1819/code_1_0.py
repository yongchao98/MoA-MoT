import sympy

def solve_pyramid_flux():
    """
    Calculates the flux of a vector field through specified sides of a pyramid.
    """
    x, y, z = sympy.symbols('x y z')
    
    # The vector field F represents the energy flow.
    F_vector_field = sympy.Matrix([3*x**3*y**2*z, 3*x**2*y**3, z])
    
    # The problem describes a square pyramid with its base on the xy-plane.
    # The sides are painted yellow and blue, interspersed. This leads to an
    # ambiguity. We will assume the "front" and "back" faces are yellow.

    # --- Calculation for Face S1 (front face, z = 4 - 4x) ---
    # The outward normal vector for this face is (4, 0, 1).
    F_on_S1 = F_vector_field.subs(z, 4 - 4*x)
    normal_S1 = sympy.Matrix([4, 0, 1])
    integrand_S1 = F_on_S1.dot(normal_S1)
    
    # The projection of S1 onto the xy-plane is a triangle defined by
    # x from 0 to 1, and y from -x to x.
    flux_S1 = sympy.integrate(integrand_S1, (y, -x, x), (x, 0, 1))

    # --- Calculation for Face S3 (back face, z = 4 + 4x) ---
    # The outward normal vector for this face is (-4, 0, 1).
    F_on_S3 = F_vector_field.subs(z, 4 + 4*x)
    normal_S3 = sympy.Matrix([-4, 0, 1])
    integrand_S3 = F_on_S3.dot(normal_S3)
    
    # The projection of S3 onto the xy-plane is a triangle defined by
    # x from -1 to 0, and y from x to -x.
    flux_S3 = sympy.integrate(integrand_S3, (y, x, -x), (x, -1, 0))
    
    # The total flux is the sum of the fluxes through the two yellow sides.
    total_flux = flux_S1 + flux_S3
    
    print("The total energy flow through the yellow sides is the sum of the individual fluxes.")
    print(f"Flux through the front face (S1) = {flux_S1}")
    print(f"Flux through the back face (S3) = {flux_S3}")
    print(f"Total Flux = {flux_S1} + {flux_S3} = {total_flux}")

solve_pyramid_flux()
<<<80/21>>>