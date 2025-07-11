import sympy

def solve_flux_integral():
    """
    Calculates the flux of the vector field F through the two yellow sides of a pyramid.
    """
    # Define the symbolic variables
    y, z = sympy.symbols('y z')

    # The problem is solved by calculating the flux through one yellow face and multiplying by 2 due to symmetry.
    # The flux integral for one face is split into two parts.
    
    # Part 1 of the integrand: 3 * x**3 * y**2 * z
    # On the surface of the face, x = 1 - z/4
    integrand1 = 3 * (1 - z/4)**3 * y**2 * z

    # Part 2 of the integrand: z/4
    integrand2 = z / 4

    # Define integration limits for the projection on the yz-plane
    # The projection is a triangle. For a given z, y ranges from -(1-z/4) to (1-z/4).
    y_limits = (y, -(1 - z/4), 1 - z/4)
    z_limits = (z, 0, 4)

    # Calculate the integral for each part
    flux_part1 = sympy.integrate(integrand1, y_limits, z_limits)
    flux_part2 = sympy.integrate(integrand2, y_limits, z_limits)

    # Calculate the total flux for one face
    flux_one_face = flux_part1 + flux_part2
    
    # The total flux is for two identical yellow faces
    total_flux = 2 * flux_one_face
    
    # Print the step-by-step calculation
    print("Based on symmetry, the total energy flow is twice the flow through one yellow face.")
    print("The flux integral for one face is calculated in two parts:")
    print(f"Integral Part 1 = {flux_part1}")
    print(f"Integral Part 2 = {flux_part2}")
    print(f"\nFlux through one yellow face = (Part 1) + (Part 2)")
    print(f"Flux through one yellow face = {flux_part1} + {flux_part2} = {flux_one_face}")
    print(f"\nTotal flux through the two yellow faces = 2 * (Flux through one face)")
    print(f"Total flux = 2 * ({flux_one_face}) = {total_flux}")

solve_flux_integral()
print(f"\n<<<80/21>>>")
