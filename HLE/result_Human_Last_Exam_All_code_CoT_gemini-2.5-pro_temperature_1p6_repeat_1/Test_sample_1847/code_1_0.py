def solve_shock_profile():
    """
    This function derives and prints the analytical solution for the density profile
    of a Mach sqrt(3) shock wave under the specified conditions.
    """

    # 1. Define the components of the final equation.
    # The density is normalized by the ambient (upstream) density rho_0.
    density_ratio = "rho/rho_0"
    
    # The position x is normalized by the characteristic length scale L.
    position_ratio = "x/L"

    # 2. The derived implicit analytical solution has the form:
    # A * (rho/rho_0 - B)^2 / (rho/rho_0 * (C - rho/rho_0)) = exp(D/E * x/L)
    # The derivation yields the following numerical coefficients:
    A = 8
    B = 1  # The upstream normalized density
    C = 2  # The downstream normalized density
    D = 4
    E = 3

    # 3. Construct the equation string, clearly showing the numbers.
    equation = (
        f"({A} * ({density_ratio} - {B})**2) / "
        f"({density_ratio} * ({C} - {density_ratio})) "
        f"= exp(({D}/{E}) * {position_ratio})"
    )

    # 4. Print the explanation and the final result.
    print("The analytical solution for the density profile is an implicit equation.")
    print("Let rho be the density and rho_0 be the ambient density.")
    print("Let x be the position and L be the ambient conductive length scale.")
    print("\nThe equation relating the dimensionless density (rho/rho_0) to the")
    print("dimensionless position (x/L) is:")
    print("\n" + equation + "\n")
    print(f"where rho/rho_0 is the density in units of the ambient density,")
    print(f"and x/L is the position in units of the ambient conductive length scale.")
    print("The numerical coefficients in the equation are:")
    print(f"Upstream normalized density = {B}")
    print(f"Downstream normalized density = {C}")
    print(f"Coefficients in the exponential term = {D}/{E}")
    print(f"Coefficient on the left side = {A}")
    
# Execute the function to get the solution
solve_shock_profile()
