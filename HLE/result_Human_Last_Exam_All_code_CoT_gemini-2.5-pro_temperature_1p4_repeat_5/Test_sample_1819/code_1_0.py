import sympy
from sympy import Symbol, integrate, Rational

def solve_flux():
    """
    This function calculates the energy flow through the yellow sides of the pyramid
    by finding the average flux and then summing it for two sides.
    """
    # Define symbols
    x, y, z = sympy.symbols('x y z')

    # Define the vector field F
    Fx = 3 * x**3 * y**2 * z
    Fy = 3 * x**2 * y**3
    Fz = z
    F = sympy.Matrix([Fx, Fy, Fz])

    # --- Step 1: Calculate flux for the front/back faces (+/- y direction) ---

    # For the front face, the plane is y = 1 - z/4.
    # The outward normal vector from projecting onto the xz-plane is N_front = (0, 1, 1/4).
    N_front = sympy.Matrix([0, 1, Rational(1, 4)])
    y_on_surface = 1 - z/4
    
    # The integrand is F · N_front, with y substituted.
    integrand_front = F.subs(y, y_on_surface).dot(N_front)

    # The projection on the xz-plane is a triangle. We define the integration limits.
    # z goes from 0 to 4. For each z, x goes from -(1-z/4) to (1-z/4).
    limit_x_lower = -(1 - z/4)
    limit_x_upper = 1 - z/4

    # Integrate to find the flux for one front face.
    flux_front_face = integrate(integrand_front, (x, limit_x_lower, limit_x_upper), (z, 0, 4))
    
    # The flux through the opposite back face is identical.
    flux_y_pair = 2 * flux_front_face

    # --- Step 2: Calculate flux for the left/right faces (+/- x direction) ---

    # For the right face, the plane is x = 1 - z/4.
    # The outward normal vector from projecting onto the yz-plane is N_right = (1, 0, 1/4).
    N_right = sympy.Matrix([1, 0, Rational(1, 4)])
    x_on_surface = 1 - z/4

    # The integrand is F · N_right, with x substituted.
    integrand_right = F.subs(x, x_on_surface).dot(N_right)

    # The projection on the yz-plane is an identical triangle. Limits for y are the same as for x.
    limit_y_lower = -(1 - z/4)
    limit_y_upper = 1 - z/4
    
    # Integrate to find the flux for one right face.
    flux_right_face = integrate(integrand_right, (y, limit_y_lower, limit_y_upper), (z, 0, 4))

    # The flux through the opposite left face is identical.
    flux_x_pair = 2 * flux_right_face

    # --- Step 3: Calculate the average and final result ---

    # Total flux through all four sides
    total_side_flux = flux_y_pair + flux_x_pair
    
    # Average flux through a single side
    avg_flux_per_side = total_side_flux / 4
    
    # The total flow through two yellow sides is twice the average per side.
    final_answer = 2 * avg_flux_per_side

    # Print the explanation and final equation
    print("Due to the asymmetry in the vector field F, the flux differs for each pair of opposite faces.")
    print(f"Flux through front/back faces: {flux_y_pair}")
    print(f"Flux through left/right faces: {flux_x_pair}")
    print("\nTo resolve the ambiguity of which pair is 'yellow', we find the average flow.")
    print(f"The average energy flow through a single side is {avg_flux_per_side}.")
    print("\nThe total energy flow through the two yellow sides is the sum of the average flow for each side:")
    print(f"{avg_flux_per_side} + {avg_flux_per_side} = {final_answer}")

solve_flux()
<<<92/21>>>