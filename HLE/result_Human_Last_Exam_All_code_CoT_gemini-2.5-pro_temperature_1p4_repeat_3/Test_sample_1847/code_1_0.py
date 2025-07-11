def solve_shock_profile():
    """
    This script presents the analytical solution for the density profile of a Mach sqrt(3) shock wave.

    The solution is derived from the Navier-Stokes equations for a monatomic ideal gas
    with a Prandtl number of 3/4. The origin of the x-coordinate is chosen at the
    point where the normalized density rho/rho_0 = 4/3.
    """

    # Define the components of the equation for clarity
    normalized_density = "rho(x) / rho_0"
    normalized_position = "(x / L)"
    exponential_term = f"exp( (4/3) * {normalized_position} )"
    
    # Construct the full equation string
    equation = f"{normalized_density} = 1 + sqrt( {exponential_term} / (8 + {exponential_term}) )"

    print("The analytical solution for the density profile is given by the following equation:")
    print("")
    print(equation)
    print("")
    print("Where:")
    print("  - rho(x) is the density at position x.")
    print("  - rho_0 is the ambient density far upstream of the shock.")
    print("  - x is the position along the direction of shock propagation.")
    print("  - L is the characteristic conductive length scale, defined as:")
    print("    L = kappa / (rho_0 * M * c_0 * C_v)")
    print("    where kappa is the thermal conductivity, M is the Mach number (sqrt(3)),")
    print("    c_0 is the ambient sound speed, and C_v is the specific heat.")
    print("  - The origin x=0 is set at the point where rho(0)/rho_0 = 4/3.")
    print("")
    print("The distinct numbers appearing in the final equation are:")
    # Using a set to show only the unique numbers
    numbers = {1, 4, 3, 8}
    for num in sorted(list(numbers)):
        print(f"  {num}")

if __name__ == '__main__':
    solve_shock_profile()
