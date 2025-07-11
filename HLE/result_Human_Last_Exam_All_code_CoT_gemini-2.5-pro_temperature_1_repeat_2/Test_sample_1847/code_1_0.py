def solve_shock_profile_equation():
    """
    This function provides the analytical solution for the density profile of a specific shock wave.

    The problem considers a Mach sqrt(3) shock in a monatomic ideal gas (gamma=5/3)
    with a Prandtl number of 3/4. This specific case, known as Becker's solution,
    is analytically solvable.

    The final result is an implicit equation relating the non-dimensional density
    to the non-dimensional position.
    """

    # Define the non-dimensional variables for clarity in the output.
    density_ratio = "rho/rho_0"
    position_ratio = "x/L"

    # The numbers in the final equation, derived from the physics.
    num1 = 1
    num2 = 2
    num3 = 1
    num4 = 3
    num5 = 8
    num6 = 3

    # Construct the equation string.
    # The origin x=0 is defined where rho/rho_0 = 1.5.
    equation = (
        f"The analytical solution is an implicit equation for the density profile:\n\n"
        f"Let rho be the density and rho_0 be the ambient density.\n"
        f"Let x be the position and L be the ambient conductive length scale L = kappa/(rho_0 * M * c_0 * C_v).\n\n"
        f"The equation is:\n"
        f"({density_ratio} - {num1})^2 / (({density_ratio}) * ({num2} - {density_ratio})) = "
        f"({num3}/{num4}) * exp(({num5}/{num6}) * {position_ratio})"
    )

    print(equation)

# Execute the function to print the solution.
solve_shock_profile_equation()
