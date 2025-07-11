import math

def print_solution():
    """
    This function prints the analytical solution for the density profile of the specified shock wave.
    The equation is derived based on the given physical parameters.
    """
    
    # The analytical solution is derived for a Mach sqrt(3) shock in a monatomic ideal gas (gamma=5/3)
    # with a Prandtl number of 3/4. The density rho is normalized by the ambient density rho_0,
    # and the position x is normalized by the characteristic length scale L, giving xi = x/L.
    
    # The derived constants in the final equation are:
    c1 = 1
    c2_num = 2
    c2_den = 3
    c3 = 8
    c4_num = 4
    c4_den = 3

    # We construct the equation string using f-strings for clarity and include Unicode characters
    # for better readability.
    
    equation_str = (
        f"ρ(ξ)/ρ₀ = {c1} + exp( ({c2_num}/{c2_den})*ξ ) / "
        f"sqrt( {c3} + exp( ({c4_num}/{c4_den})*ξ ) )"
    )
    
    explanation_str = (
        "where:\n"
        "  - ρ(ξ) is the density at the non-dimensional position ξ.\n"
        "  - ρ₀ is the ambient (pre-shock) density.\n"
        "  - ξ = x/L is the non-dimensional position.\n"
        "  - exp(z) is the exponential function e^z, and sqrt(z) is the square root of z."
    )

    print("The analytical solution for the density profile is:")
    print(equation_str)
    print("\n" + explanation_str)

# Execute the function to print the solution
print_solution()