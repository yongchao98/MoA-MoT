import sympy as sp

def solve_flux_problem():
    """
    This function calculates the total energy flux through the two yellow sides of a pyramid.
    """
    # Define symbolic variables
    x, y, z = sp.symbols('x y z')

    # Define the vector field F = (3x³y²z, 3x²y³, z)
    F = sp.Matrix([3 * x**3 * y**2 * z, 3 * x**2 * y**3, z])

    # --- Calculation for Yellow Side 1 (in the y > 0 region) ---
    # The plane containing this side is 4y + z = 4. We can write y in terms of z.
    y_S1 = (4 - z) / 4
    
    # The outward normal vector N1 is (0, 4, 1).
    N1 = sp.Matrix([0, 4, 1])
    
    # The dot product F · N1
    F_dot_N1 = F.dot(N1)
    
    # Substitute y from the plane equation into the dot product to get the integrand.
    integrand = F_dot_N1.subs(y, y_S1)
    
    # We project onto the xz-plane. The integration bounds for z are from 0 to 4.
    # For a given z, the bounds for x are from -(1-z/4) to (1-z/4).
    x_lower_bound = -(1 - z/4)
    x_upper_bound =  (1 - z/4)
    
    # First, integrate with respect to x.
    integral_over_x = sp.integrate(integrand, (x, x_lower_bound, x_upper_bound))
    
    # Now, integrate the result with respect to z from 0 to 4.
    # This gives the value of the numerator in the flux formula: ∫∫_R (F · N) dx dz
    flux_numerator_1 = sp.integrate(integral_over_x, (z, 0, 4))
    
    # The flux formula is Flux = ∫∫ (F · N) dx dz / |N · j|, where j = (0,1,0).
    # |N1 · j| = |(0,4,1) · (0,1,0)| = 4.
    denominator_1 = abs(N1.dot(sp.Matrix([0, 1, 0])))
    
    flux_1 = flux_numerator_1 / denominator_1

    # --- Calculation for Yellow Side 2 (in the y < 0 region) ---
    # Due to the symmetries in the pyramid's geometry and the vector field components,
    # the flux through the second yellow side is identical to the first.
    flux_2 = flux_1

    # --- Total Flux ---
    total_flux = flux_1 + flux_2

    # --- Print the results step-by-step ---
    print("The energy flow (flux) is the surface integral of the vector field F over each yellow side.")
    print("\nCalculating the flux for the first yellow side (y > 0):")
    print(f"The integral ∫∫ (F · N1) dA evaluates to {flux_numerator_1}.")
    print(f"The normalization factor |N1 · j| is {denominator_1}.")
    print(f"Flux_1 = {flux_numerator_1} / {denominator_1} = {flux_1}")

    print("\nThe flux for the second yellow side (y < 0) is identical due to symmetry.")
    print(f"Flux_2 = {flux_2}")
    
    print("\nThe total energy flow through both yellow sides is the sum of their individual fluxes:")
    print(f"Total Flux = {flux_1} + {flux_2} = {total_flux}")
    
    # The numerical value is approx 4.952
    print(f"\nNumerical Answer: {total_flux.evalf()}")


solve_flux_problem()
<<<104/21>>>